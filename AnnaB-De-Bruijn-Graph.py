#De Bruijn Graph Simulator


import random
#import networkx as nx
#import matplotlib.pyplot as plt

def randomDNA(length,nucleotide_frequencies):
#Generate a random DNA text
    nucleotides= ['a', 'c', 't', 'g']
    weights= [i/sum(nucleotide_frequencies) for i in nucleotide_frequencies]
    random_sequence=random.choices(nucleotides, weights=weights, k=length)
    random_sequence= "".join(random_sequence)

    return random_sequence




def add_repetitive_sequences(dna_sequence, repeat_length, num_repeats):
#Add repetitive elements
    repetitive_sequence=''
    for i in range(repeat_length):
        var=random.choice('actg')
        repetitive_sequence = repetitive_sequence + var
    for i in range(num_repeats):
        insertion_point = random.randint(0, len(dna_sequence))
        dna_sequence = dna_sequence[:insertion_point] + repetitive_sequence + dna_sequence[insertion_point:]
    return dna_sequence



def simulate_sequencing(dna_sequence, read_length, coverage, error_rate):
    reads = []
    for i in range(int(coverage * len(dna_sequence) / read_length)):
        start_position = random.randint(0, len(dna_sequence) - read_length)
        read = dna_sequence[start_position:start_position + read_length]

        # Introduce errors based on the error rate
        for i in range(read_length):
            if random.random() < error_rate:
                read = read[:i] + random.choice('actg'.replace(read[i], '')) + read[i+1:]
        reads.append(read)
    return reads


def kmers(kvalue,sequencing_reads):
#From the reads, generate kmers 
    kmerlist=[]
    for i in range(0,len(sequencing_reads)):
        read=sequencing_reads[i]
        for j in range(0,int(len(read)/kvalue)):
            kmer1=read[j:j+kvalue]
            kmerlist.append(kmer1)
        
    return kmerlist

def construct_de_bruijn_graph(kmer_list):
    graph = {}
    for kmer in kmer_list:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if prefix not in graph:
            graph[prefix] = [suffix]
        else:
            if suffix not in graph[prefix]:
                graph[prefix].append(suffix)

    return graph

def find_contig_paths(graph, node, visited, current_path):
    visited[node] = True
    all_contig_paths = []
    current_path.append(node)

    if node in graph:
        for neighbor in graph[node]:
            #print('neighbor', neighbor)
            if neighbor not in visited or visited[neighbor]: #this is the last one, either the node does not exist in the graph or it has been visited
                # here we don't want to add this last node to current_path because it will effect the other iterations on the loop
                final_path = current_path[:]
                final_path.append(neighbor)
                all_contig_paths.append(final_path)
            else:
                new_visited = dict(visited) # we want to make this variable specific to the path we are building so we create a new one
                contig_paths = find_contig_paths(graph, neighbor, new_visited, current_path)
                all_contig_paths += contig_paths
    return all_contig_paths
               
def trace_contigs(graph):
    visited = {node: False for node in graph}
    all_contigs = []

    for node in graph:
        if not visited[node]:
            current_path = []
            contig_paths = find_contig_paths(graph, node, visited, current_path)
            for contig_path in contig_paths:
                # Start with the first node and then concatenate the last character of the rest of the nodes
                contig = contig_path[0]
                for node in contig_path[1:]:
                    contig += node[-1]
                all_contigs.append(contig)

    return all_contigs



def main():
    #Example variables
    nucleotide_frequencies= [0.25,0.25,0.25,0.25]
    # percents of a,c,t,g respectively
    sequence_length = 50
    random_dna = randomDNA(sequence_length,nucleotide_frequencies)

    # Bonus: Simulate addition of repetitive sequences
    repeat_length = 5
    num_repeats = 3
    dna_with_repeats = add_repetitive_sequences(random_dna, repeat_length, num_repeats)

    #Simulate a sequencing project
    read_length = 10
    coverage = 5
    error_rate = 0.02
    sequencing_reads = simulate_sequencing(dna_with_repeats, read_length, coverage, error_rate)

    #Create kmers
    kvalue=5
    kmerlist=kmers(kvalue,sequencing_reads)
    print("kmerlist:",kmerlist)

    #Create the graph 
    de_bruijn_graph = construct_de_bruijn_graph(kmerlist)

    print('De Bruijn graph')
    for node, neighbors in de_bruijn_graph.items():
        print(f"{node} -> {','.join(neighbors)}")
    print('--------------------------------------------')
    #print('Walking the De Bruijn graph')
    contigs = trace_contigs(de_bruijn_graph)

    #Print the resulting contigs
    for i, contig in enumerate(contigs, 1):
        print(f"Contig {i}: {contig}")
        
    #draw_de_bruijn_graph(de_bruijn_graph) 
    
    
main()


     
