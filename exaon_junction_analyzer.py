import pysam
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import sys
import os

try:
    import torch
    from torch_geometric.data import Data
    TORCH_AVAILABLE = True
    print("true")
except ImportError:
    TORCH_AVAILABLE = False
    print("false")

class SpliceGraphAnalyzer:
    def __init__(self, bam_file, gff_file):
        self.bam_file = bam_file
        self.gff_file = gff_file
        
        #junction counts
        self.junctions = Counter()
        
        #exon intervals
        self.known_exons = defaultdict(list)
        
        self.graph = nx.DiGraph()
        self.chrom_map = {} 

    def load_annotation(self):
        count = 0
        try:
            with open(self.gff_file, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    parts = line.strip().split('\t')
                    if len(parts) < 9: continue
                    
                    feature_type = parts[2]
                    chrom = parts[0]
                    
                    if not chrom.startswith('chr'):
                        chrom = 'chr' + chrom
                    
                    if chrom not in self.chrom_map:
                        self.chrom_map[chrom] = len(self.chrom_map)

                    if feature_type == 'exon':
                        start = int(parts[3])
                        end = int(parts[4])
                        
                        self.known_exons[chrom].append((start - 1, end))
                        count += 1
                        
            print(f"{count} exons from reference.")
            
        except Exception as e:
            print(f"Error reading annotation file: {e}")
            sys.exit(1)

    def process_bam(self):
        #does NOT require a .bai index file
        try:
            samfile = pysam.AlignmentFile(self.bam_file, "rb", check_sq=False)
            
            read_count = 0
            
            for read in samfile:
                read_count += 1
                if read.is_unmapped or read.is_secondary:
                    continue

                if not read.cigarstring or 'N' not in read.cigarstring:
                    continue

                chrom = read.reference_name
                
                if chrom not in self.chrom_map:
                     self.chrom_map[chrom] = len(self.chrom_map)

                current_pos = read.reference_start
                
                if read.cigartuples:
                    for op, length in read.cigartuples:
                        if op in [0, 2, 7, 8]: 
                            current_pos += length
                        
                        elif op == 3: 
                            j_start = current_pos
                            j_end = current_pos + length
                          
                            self.junctions[(chrom, j_start, j_end)] += 1
                            
                            current_pos += length
                        
            samfile.close()
            print(f"Read counts {read_count} reads.")
            print(f"{len(self.junctions)} unique splice junctions.")
            
        except Exception as e:
            print(f"Error reading BAM file: {e}")
            import traceback
            traceback.print_exc()

    def build_graph(self):
        #Nodes: Genomic coordinates (e.g., 'chr1:100')
        #Edges: Splicing events       
        for (chrom, start, end), count in self.junctions.items():
            u = f"{chrom}:{start}" # Donor Node
            v = f"{chrom}:{end}"   # Acceptor Node
            
            self.graph.add_edge(u, v, weight=count)
            
            is_known = self._check_is_known(chrom, start, end)
            
            self.graph[u][v].update({
                'type': 'known' if is_known else 'novel',
                'chrom': chrom,
                'start': start,
                'end': end
            })

    def _check_is_known(self, chrom, j_start, j_end):
        if chrom not in self.known_exons: 
            return False
        
        start_match = False
        for _, e_end in self.known_exons[chrom]:
            if e_end == j_start:
                start_match = True
                break
        
        end_match = False
        if start_match:
            for e_start, _ in self.known_exons[chrom]:
                if e_start == j_end:
                    end_match = True
                    break
        
        return start_match and end_match

    def find_lsvs(self, min_reads=5, min_psi=0.10):
        #Local Splice Variations (LSVs).

        print(f"Min Reads={min_reads}, Min PSI={min_psi*100}%)...")
        lsvs = []
        
        for node in self.graph.nodes():
            out_edges = list(self.graph.out_edges(node, data=True))
            in_edges = list(self.graph.in_edges(node, data=True))
            
            if len(out_edges) > 1:
                total_reads = sum([d['weight'] for u, v, d in out_edges])
                
                if total_reads >= min_reads:
                    lsv_data = {
                        'type': 'Source_LSV',
                        'ref_node': node,
                        'total_reads': total_reads,
                        'edges': []
                    }
                    
                    valid_lsv = False
                    for u, v, d in out_edges:
                        psi = d['weight'] / total_reads
                        if psi >= min_psi: valid_lsv = True
                        
                        lsv_data['edges'].append({
                            'target': v, 
                            'reads': d['weight'], 
                            'psi': psi, 
                            'annotation': d['type'],
                            'coords': (d['chrom'], d['start'], d['end'])
                        })
                    
                    if valid_lsv: lsvs.append(lsv_data)

            if len(in_edges) > 1:
                total_reads = sum([d['weight'] for u, v, d in in_edges])
                
                if total_reads >= min_reads:
                    lsv_data = {
                        'type': 'Target_LSV',
                        'ref_node': node,
                        'total_reads': total_reads,
                        'edges': []
                    }
                    
                    valid_lsv = False
                    for u, v, d in in_edges:
                        psi = d['weight'] / total_reads
                        if psi >= min_psi: valid_lsv = True
                        
                        lsv_data['edges'].append({
                            'source': u, 
                            'reads': d['weight'], 
                            'psi': psi, 
                            'annotation': d['type'],
                            'coords': (d['chrom'], d['start'], d['end'])
                        })
                    
                    if valid_lsv: lsvs.append(lsv_data)
                    
        lsvs.sort(key=lambda x: x['total_reads'], reverse=True)
        return lsvs

    def visualize_lsv(self, lsv_data, lsv_id):
        subgraph = nx.DiGraph()
        labels = {}
        colors = []
        pos = {}
        
        ref_node = lsv_data['ref_node']
        ref_pos = int(ref_node.split(':')[1])
        
        subgraph.add_node(ref_node)
        pos[ref_node] = (ref_pos, 0)
        
        num_edges = len(lsv_data['edges'])
        if num_edges > 1:
            y_offsets = [i - (num_edges - 1) / 2 for i in range(num_edges)]
            y_offsets.reverse() 
        else:
            y_offsets = [0]

        for i, edge in enumerate(lsv_data['edges']):
            other_node = edge.get('target', edge.get('source'))
            other_pos_val = int(other_node.split(':')[1])
            
            subgraph.add_node(other_node)
            
            pos[other_node] = (other_pos_val, y_offsets[i] * 50)
            
            if lsv_data['type'] == 'Source_LSV':
                subgraph.add_edge(ref_node, other_node)
            else:
                subgraph.add_edge(other_node, ref_node)
            
            key = (ref_node, other_node) if lsv_data['type'] == 'Source_LSV' else (other_node, ref_node)
            labels[key] = f"{edge['psi']*100:.1f}%"
            colors.append('blue' if edge['annotation'] == 'known' else 'red')

        plt.figure(figsize=(12, 6))
        plt.title(f"LSV #{lsv_id} ({lsv_data['type']}) - Reads: {lsv_data['total_reads']}")
        
        nx.draw(subgraph, pos, with_labels=True, node_color='lightgray', 
                edge_color=colors, width=2.5, node_size=3500, font_size=9, 
                font_weight='bold', connectionstyle='arc3,rad=0.1')
        
        nx.draw_networkx_edge_labels(subgraph, pos, edge_labels=labels, font_color='black', label_pos=0.5)
        
        plt.margins(0.2)
        
        filename = f"lsv_{lsv_id}.png"
        plt.savefig(filename, dpi=150)
        plt.close()

    def export_to_bed(self, lsv_results, filename="lsv_results.bed"):
        try:
            with open(filename, 'w') as f:
                f.write(f"track name='LSV_Analysis' description='Splice Junctions' itemRgb='On'\n")
                
                for i, lsv in enumerate(lsv_results):
                    for edge in lsv['edges']:
                        chrom, start, end = edge['coords']
                        
                        color = "0,0,255" if edge['annotation'] == 'known' else "255,0,0"
                        
                        name = f"LSV{i+1}|psi:{edge['psi']:.2f}"
                        score = int(edge['psi'] * 1000)
                        
                        line = f"{chrom}\t{start}\t{end}\t{name}\t{score}\t.\t{start}\t{end}\t{color}\n"
                        f.write(line)
        except Exception as e:
            print(f"Error writing BED file: {e}")

    def export_to_pyg(self, output_file="splice_graph.pt"):
        if not TORCH_AVAILABLE:
            print("Please install the pytorch libraries.")
            return

        nodes = list(self.graph.nodes())
        node_to_idx = {node: i for i, node in enumerate(nodes)}
        
        node_features = []
        for node in nodes:
            chrom_str, pos_str = node.split(':')
            chrom_id = self.chrom_map.get(chrom_str, -1)
            pos_val = float(pos_str)
            node_features.append([chrom_id, pos_val])
            
        x = torch.tensor(node_features, dtype=torch.float)

        edge_indices = []
        edge_attrs = []
        
        for u, v, data in self.graph.edges(data=True):
            src = node_to_idx[u]
            dst = node_to_idx[v]
            edge_indices.append([src, dst])
            
            weight = float(data['weight'])
            is_known = 1.0 if data['type'] == 'known' else 0.0
            edge_attrs.append([weight, is_known])

        if len(edge_indices) == 0:
            print("Graph is empty. Nothing to save.")
            return

        edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_attrs, dtype=torch.float)

        data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
        
        data.num_nodes = len(nodes)
        
        torch.save(data, output_file)
        print(f"output file: {output_file}")
        print(f"{data.num_nodes} Nodes, {data.num_edges} Edges.")
        print(f"x={data.x.shape}, edge_attr={data.edge_attr.shape}")

if __name__ == "__main__":
    BAM_PATH = "short_reads.bam"
    GFF_PATH = "human-chr19_P.gff"
    
    MIN_READS = 5
    MIN_PSI = 0.10
    
    if os.path.exists(BAM_PATH) and os.path.exists(GFF_PATH):
        analyzer = SpliceGraphAnalyzer(BAM_PATH, GFF_PATH)
        analyzer.load_annotation()
        analyzer.process_bam()
        analyzer.build_graph()
        lsv_results = analyzer.find_lsvs(min_reads=MIN_READS, min_psi=MIN_PSI)
        
        print(f"{len(lsv_results)} significant LSVs.")
        for i, lsv in enumerate(lsv_results[:5]):
            print(f"LSV #{i+1} | Type: {lsv['type']}")
            print(f"Reference Node: {lsv['ref_node']} | Total Support: {lsv['total_reads']} reads")
            
            for edge in lsv['edges']:
                target = edge.get('target', edge.get('source'))
                print(f"Connects to: {target}")
                print(f"PSI: {edge['psi']:.2f} ({edge['psi']*100:.1f}%) | Annotation: {edge['annotation']}")
            
            analyzer.visualize_lsv(lsv, i+1)

        if lsv_results:
            analyzer.export_to_bed(lsv_results, filename="lsv_results.bed")
        analyzer.export_to_pyg(output_file="splice_graph.pt")    
 
