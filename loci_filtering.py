#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import sys
import argparse
from homolog4 import *
from collections import defaultdict

# Copyright(C) 2014 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description="Filter out redundant hits, by loci, from the initial BLAST parse and remove organisms that lack neighborhoods.")
                
    parser.add_argument("-i", "--infolder", dest="infolder", default='./blast_parse_raw_operon/', metavar="FOLDER",
                help="A folder that contains the initial parse of the BLAST hits. This program assumes that no loci filtering or organism removal has been done yet.")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="FOLDER", default='./blast_parse/',
                help="Folder where the BLAST results will be stored. Default is the folder './blast_result/'.")
    
    parser.add_argument("-q", "--operon_query", dest="operon_query", default='./regulonDB/operon_names_and_genes.txt', metavar="FILE",
                help="A file that contains the names and genes comprising the operons that are under investigation.")

    parser.add_argument("-r", "--reference", dest="reference", default='NC_000913', metavar="STRING",
                help="An accession number of the reference organism.")
    
    parser.add_argument("-f", "--filter", dest="filter", default='', metavar="FILE",
                help="A file that contains the accession numbers of the organisms that are under investigation.")            
    
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
    
    parser.add_argument("-g", "--max_gap", dest="max_intergenic_gap", metavar="INT", default = 500, type=int,
                help="Length of the largest allowable intergenic gap allowed in determining a gene neighborhood. Default is 500 nucleotides")
    
                
    return parser.parse_args()


def check_options(parsed_args):
    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder
    else:
        print "The folder %s does not exist." % parsed_args.infolder
        sys.exit()
    
    # if the directory that the user specifies does not exist, then the program makes it for them. 
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    outfolder = parsed_args.outfolder
    if outfolder[-1] != '/':
        outfolder = outfolder + '/'
        
    if os.path.exists(parsed_args.operon_query):
        operon_query = parsed_args.operon_query
    else:
        print "The file %s does not exist." % parsed_args.operon_query
        sys.exit()

    if os.path.exists(parsed_args.filter):
        filter_file = parsed_args.filter
    elif parsed_args.filter == '':
        filter_file = parsed_args.filter
    else:
        print "The file %s does not exist." % parsed_args.filter
        sys.exit()
        
    if os.path.exists(parsed_args.operon_query):
        filter_file = parsed_args.operon_query
    else:
        print "The file %s does not exist." % parsed_args.operon_query
        sys.exit()
        
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)
    
    
    if parsed_args.num_proc < 1:
        max_intergenic_gap = 1
    else:
        max_intergenic_gap = int(parsed_args.max_intergenic_gap)    

    return infolder, outfolder, operon_query, filter_file, num_proc, operon_query, max_intergenic_gap

#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result

# this function will return a dictionary of operon keyed off the operon name with data values in the form
# of a list of homologs which are homologous. ex. [abcA, abcB]
def return_self_homolog_dict(operon_list = 'operon_name_and_genes.txt', prot_file = 'operon_protein_query.fa', rna_file = 'operon_rna_query.fa'):
    # makes a dictionary keyed by operon name and a list of the gene contained by the operon
    operon_dict = {}
    for line in [i.strip() for i in open(operon_list).readlines()]:
        tmp = line.split('\t')
        operon_dict.update({tmp[0]:tmp[1:]})
        
    # set up databases for the different types of genes
    # for proteins -p must be set to true    
    cmd = "formatdb -i %s -p T -o F" % (prot_file)
    os.system(cmd)
    # for RNA genes -p must be set to false 
    #cmd = "formatdb -i %s -p F -o F" % (rna_file)
    #os.system(cmd)
    
    # blast each set of genes against itself
    '''cmd = "blastall -p blastp -a %i -i %s -d %s -e %s -o %s -m 9" % (os.sysconf("SC_NPROCESSORS_ONLN"), prot_file, prot_file, '1e-10', 'self_prot.txt')
    os.system( cmd )'''
    cmd = "blastall -p blastp -a %i -i %s -d %s -e %s -o %s -m 8" % (os.sysconf("SC_NPROCESSORS_ONLN"), prot_file, prot_file, '1e-10', 'self_prot.txt')
    os.system( cmd )
    #cmd = "blastall -p blastn -a %i -i %s -d %s -e %s -o %s -m 9" % (os.sysconf("SC_NPROCESSORS_ONLN"), rna_file, rna_file, '1e-10', 'self_rna.txt')
    #os.system( cmd )
    
    # in this next section i will read in the resulting blast results, and construct a dictionary which will be keyed off gene name and provide a list
    # of homologs from the operon set. This list will allow the program to filter out spurious results. We will miss fusions of homologous genes, but
    # hopefully this will be a rare event in our dataset, untill this can be revised
    
    lst = [i.strip() for i in open('self_prot.txt').readlines() if i[0] != '#']
    #for line in [i.strip() for i in open('self_rna.txt').readlines() if i[0] != '#']:
    #    lst.append(line)
        
    result = {} 
    
    print "got here 1"
    for line in lst:
        source, hit = line.split('\t')[0:2]
        source_annotation = source.split('|')[2]
        hit_annotation = hit.split('|')[2]
        # we have two genes in the test set that are homologous
        if source_annotation != hit_annotation:
            if source_annotation not in result.keys():
                result.update({source_annotation: [hit_annotation]})
            else:
                result[source_annotation].append(hit_annotation)
    print "got here 2"    
    return result
    
    
# The purpose of this function is to filter out the spurious hits on a locus, and determine the annotation of the gene at 
# that position.  To do this the program will make a dict of each [locus/start] ? and then determine the annotations that
# exist for it.  If there are two annotations that have homologous genes then the best hit will be used. if there are two
# annotations for a locus which are not homologous then some sort of hit analysis will be performed to determine is there
# is a good candidate for a gene fusion.  (should look at papers on this).  When done the function will report a list of 
# homologs that are ordered by start position whichi have been filtered for the best hit, or as a fusion.
def filter_locus_hits(h_list, self_homolog_dict):
    pass    
    
def main():
    start = time.time()
    
    print time.time() - start

    # ./blast_parse.py -f phylo_order.txt
if __name__ == '__main__':
    main()    
    
    
    
    
    
    
    
    
    
    
    
    

