#!/usr/bin/python

import time
import os
import sys
import argparse

#TODO:  Rename this script, it's horrible!

# Copyright(C) 2014 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():
    
    parser = argparse.ArgumentParser(description='The purpose of this script is to run the full software suite that we have developed to study operons using as few inputs as possible.  This will facilitate the ease of use as much as possible.')

    parser.add_argument("-i", "--infile", dest="infile", metavar="FILE", default='',
                help="FILE containing the input for the pipeline stage chosen. Most stages require a FOLDER containing the files the pipeline will work on. The option -I should be used.")
                
    parser.add_argument("-I", "--infolder", dest="infolder", metavar="FOLDER", default='/home/dave/Desktop/all_genbank',
                help="Folder containing all files used by the pipeline stage. This is the most common way to start the program.")
                 
    parser.add_argument("-o", "--outfile", dest="outfile", metavar="FILE", default='./operon_query.fa',
                help="Folder where the final result of operon analysis is stored.")
                
    parser.add_argument("-f", "--filter", dest="filter", metavar="FILE", default='',
                help="File restrictiong which accession numbers this script will process. If no file is provided, filtering is not performed.")
                
    parser.add_argument("-p", "--operon_file", dest="operon_file", metavar="FILE", default='./regulonDB/operon_names_and_genes.txt',
                help="File which contains operon information for use in operon queries. The file format is operon_name followed by the constituent gene names, tab delineated.")
                
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")

    parser.add_argument("-r", "--refrence", dest="refrence", metavar="STRING", default = 'NC_000913',
                help="Accession number of the refrence organism. This information is used to determine the product type of each gene (RNA/Protein), a necessary piece of information to classify the operons that are under investigation.")
    
    parser.add_argument("-R", "--refrence_file", dest="refrence_file", metavar="FILE", default = '',
                help="File containing a list of accession number(s) of refrence organism(s) . This information is used to generate a list of query sequences for each operon gene using CD-HIT and some intelligent filtering to disreguard possible misannotations.")
                
    parser.add_argument("-s", "--stage", dest="stage", metavar="INT", default=1,
                help="The stage in the pipeline where you wish to start from. Ususally you want everything to run from scratch, and the default will be fine. Fill in the rest later.")
                       
    return parser.parse_args()
    
    
def check_options(parsed_args):
    # section of code that checks the infolder entry    
    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder
    else:
        print "The folder %s does not exist." % parsed_args.infolder
        sys.exit()
        
    if infolder[-1] != '/':
        infolder = infolder + '/'
        
    # I'm not issuing a message that hey, this file is there and will be overwritten, the program will just overwrite.
    # TODO: check that the folder that this is saved in exists
    outfile = parsed_args.outfile
    
    # Check the operon file
    if os.path.exists(parsed_args.operon_file):
        operon_file = parsed_args.operon_file
    else:
        print "The operon file %s does not exist." % parsed_args.operon_file
        sys.exit()
    
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)

    # Code block that determines what the refrence organism(s) are, and return the result as a list.
    if parsed_args.refrence_file == '' or os.path.exists(parsed_args.refrence_file):
        refrence_file = parsed_args.refrence_file
    else:
        print "The file %s does not exist." % parsed_args.refrence_file
        sys.exit()
    
    if refrence_file == '':
        refrence_list = [parsed_args.refrence]
    else:
        refrence_list = [i.strip() for i in open(parsed_args.refrence_file).readlines()]

    return infolder, outfile, operon_file, num_proc, refrence_list
    
    
def main():
    
    start = time.time()

    parsed_args = parser_code()
    
    infolder, outfile, operon_file, num_proc, refrence_list = check_options(parsed_args)
    
    print infolder, outfile, operon_file, num_proc, refrence_list
    
    print time.time() - start
    
if __name__ == '__main__':
    main()
