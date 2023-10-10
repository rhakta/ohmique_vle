#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Vincent Le"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Vincent LE"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Vincent LE"
__email__ = "le75.vincent@gmail.com"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    with open (fastq_file) as f:
        file_read=[]
        for line in f:
            file_read.append(line.strip())
        for read in file_read[1::4]:
            yield(read)
    pass


def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of size kmer_size.
    """
    for i in range(len(read)-kmer_size+1):
        k_mer=""
        for j in range(kmer_size):
            k_mer+=read[i+j]
        yield(k_mer)
    pass


def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    dictio=dict()
    for read in read_fastq(fastq_file):
        for k_mer in cut_kmer(read,kmer_size):
            try:
                dictio[k_mer]+=1
            except:
                dictio[k_mer]=1
    return(dictio)


def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    graph=nx.DiGraph()
    for key,value in kmer_dict.items():
        graph.add_node(key[0:len(key)-1])
        graph.add_node(key[1:len(key)])
        graph.add_edge(u_of_edge=key[0:len(key)-1],v_of_edge=key[1:len(key)],weight=value)
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list:
        print(path)
        graph.remove_nodes_from(path[1:-1])
        if (delete_entry_node):
            graph.remove_node(path[0])
        if delete_sink_node:
            graph.remove_node(path[-1])
    return graph


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    #check du max weight
    std_weight=statistics.stdev(data=weight_avg_list)
    if (std_weight>0):
        max_weight=max(weight_avg_list)
        num_best_path=weight_avg_list.index(max_weight)
    else:
        #check du max len
        max_length=max(path_length)
        if (max_length>0):
            num_best_path=path_length.index(max_length)
        else:
            num_best_path=randint(0,len(path_list)-1)
    del path_list[num_best_path]
    graph2=remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph2

def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    G_paths=nx.all_simple_paths(graph,source=ancestor_node,target=descendant_node)

    L_paths=[]
    for path in G_paths:
        L_paths.append(path)

    L_avg_weight=[]
    for path in L_paths:
        L_avg_weight.append(path_average_weight(graph,path))

    L_len_paths=[]
    for path in L_paths:
        L_len_paths.append(len(path))

    return (select_best_path(graph, L_paths, L_len_paths, L_avg_weight, 
                     delete_entry_node=False, delete_sink_node=False))


def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble=False
    for node in graph:
        G_predecessors=graph.predecessors(node)
        L_pred=[]
        for pred in G_predecessors:
            L_pred.append(pred)
        if (len(L_pred)>1):
            for i in range (len(L_pred)-1):
                for j in range (i+1,len(L_pred)):
                    node_ancestor=nx.lowest_common_ancestor(graph,L_pred[i],L_pred[j])
                    if (node_ancestor!=None):
                        bubble=True
                        break
        if bubble:
            break
    if bubble:
        graph=simplify_bubbles(solve_bubble(graph,node_ancestor,node))
    return graph

def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    start_nodes=[]
    for node in graph:
        if (all(False for _ in (graph.predecessors(node)))):
            start_nodes.append(node)
    return start_nodes

def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    start_nodes=[]
    for node in graph:
        if (all(False for _ in (graph.successors(node)))):
            start_nodes.append(node)
    return start_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    L_contigs=[]
    for st_node in starting_nodes:
        for end_node in ending_nodes:
            all_path=nx.all_simple_paths(graph,source=st_node,target=end_node)
            for path in all_path:
                sequ=path[0][0]
                for k_mer in path:
                    sequ+=k_mer[1]
                L_contigs.append((sequ,len(path)+1))
                print(sequ)
    return (L_contigs)

def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    num_contig=0
    str_file=""
    for contig,len_contig in contigs_list:
        str_file+=">contig_"+str(num_contig)+" len="+str(len_contig)+"\n"+contig+"\n"
        num_contig+=1
    str_file.strip()
    new_file=open(output_file,'w')
    new_file.write(str_file)
    new_file.close()
    pass


def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == '__main__': # pragma: no cover
    main()
