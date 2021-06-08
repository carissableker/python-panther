#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
import pandas as pd
import os
import sys
import argparse
import panther_api
import time

def main(inputfile, outputprefix, remove_version, organism, test_type, annotation_option, min_size, start_cluster):
    
    with open(inputfile, 'r') as p:
        # get to start line
        j = 0
        while j < start_cluster:
            p.readline()
            j += 1

        # go
        i = start_cluster
        for line in p:
            print("\nCluster %i"%i)
            genes = line.rstrip().split('\t')
            if len(genes) < min_size:
                print("Cluster too small.")
                print("-----------------------------")
                i += 1
                continue

            if remove_version:
                genes = [gene.split('.')[0] for gene in genes]
            gene_list = outputprefix + '_cluster_%i.txt'%i
            with open(gene_list, 'w') as o:
                for gene in genes:
                    o.write('%s\n'%gene)

            df = panther_api.panther_api_overrepresentation(gene_list, organism, annotation_option, test_type)
            if df is not None: 
                outputfile = outputprefix + '_cluster_%i.panther'%i
                df.to_csv(outputfile, sep='\t')
            else:
                print("No results to save")
            print("-----------------------------")
            i += 1
            time.sleep(2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PANTHER overenrichment test on the output of a clustering algorithm.')
    parser.add_argument("inputfile", type=str, help="Cluster file. One cluster per line, tab seperated gene IDs.")
    parser.add_argument("outputprefix", type=str, help="Prefix to give to output gene lists and enrichment result files. Output will be saved to OUTPUTPREFIX_cluster_<num>.txt and OUTPUTPREFIX_cluster_<num>.panther")
    parser.add_argument("--remove_version", action='store_true', default=False, help="Whether to remove version number from gene IDs")
    parser.add_argument("--organism", type=str, help="Organism for reference/background.", default='Homo sapien')
    parser.add_argument("--test_type", type=str, help="One of FISHER or BINOMIAL", default='FISHER')
    parser.add_argument("--annotation_option", type=str, help="Annotation option, see table in code", default='fullgo_bp_comp')
    parser.add_argument("--min_size", type=int, help="Minimum cluster size to test", default=3)
    parser.add_argument("--start_cluster", type=int, default=0, help="Line in file to start at (from 0)")
    args = parser.parse_args()
    
    main(args.inputfile, args.outputprefix, args.remove_version, args.organism, args.test_type, args.annotation_option, args.min_size, args.start_cluster)
