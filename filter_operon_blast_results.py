#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import sys
import argparse
import math
from homolog4 import *
from collections import defaultdict
import itertools

# Copyright(C) 2014 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description="This program will be used to remove spurious results from a BLAST search organized by operon.  The result will be the basis of the database for GBEER and for other anaylsis.")
                
    ##parser.add_argument("-i", "--infolder", dest="infolder", default='./blast_result/', metavar="FOLDER",
    ##           help="A folder that contains the operon BLAST results.")
    
    parser.add_argument("-i", "--infolder", dest="infolder", default='./intermediate_for_debug/unfiltered_operon/', metavar="FOLDER",
                help="A folder that contains the operon BLAST results.")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="FOLDER", default='./optimized_operon/',
                help="Folder where the filtered results will be stored. Default is the folder './optimized_operon/'.")

    parser.add_argument("-f", "--filter", dest="filter", default='', metavar="FILE",
                help="A file that contains the operons that are under investigation.  All others will be omitted from analysis an results.")            
    
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
                
    parser.add_argument("-e", "--eval", dest="eval", default='1e-10', metavar="FLOAT", type=float,
                help="Use this option to change the eval for the BLAST search that is permitted. Useful if you would like to investigate what altering the eval threshold will do to your results.")
                
    parser.add_argument("-g", "--max_gap", dest="max_gap", metavar="INT", default = 500, type=int,
                help="Largest allowable gap to be considered a gene block by the analysis.")
                
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
     
    if os.path.exists(parsed_args.filter):
        filter_file = parsed_args.filter
    elif parsed_args.filter == '':
        filter_file = parsed_args.filter
    else:
        print "The file %s does not exist." % parsed_args.filter
        sys.exit()
        
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)
    
    # validate the input for the eval
    try:    
        e_val = float(parsed_args.eval)
    except:
        print "The e-value you entered is not a floating point number, please enter a floating point number, ex. '1e-3', or '12'."
        sys.exit()
        
    # validate the input for the maximum allowed gap
    try:    
        max_gap = int(parsed_args.max_gap)
        if max_gap <= 0:
           print "The gap that you entered %s is a negative number, please enter a positive integer." % parsed_args.max_gap
           sys.exit()
        else:
           pass
    except:
        print "The gap that you entered %s is not an integer, please enter a positive integer." % parsed_args.max_gap
        sys.exit()
    
    return infolder, outfolder, filter_file, num_proc, e_val, max_gap


#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def return_recursive_dir_files(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# fast implementation of an order preserving make unique function
def MakeUnique(lst, function):#lambda a: a.start): 
    seen = {}
    result = []
    for item in lst:
        marker = function(item)
        if marker not in seen:
            seen.update({marker:1})
            result.append(item)
    return result


# WARNING!!!!!!! NOT, not, NOT rigorously tested yet!!!!!!! (but appears to be working just fine)
# Return a list of lists of homologs = [[],[]], number of splits, and number of duplications. 
# unique_genes_in_organism, len_operon are integers grouped_list is a list of lists, and only contains groups 2 or more.
def OptimalSplitting(grouped_lists, unique_genes_in_organism):
    len_unique_grouped = len(MakeUnique([item for sublist in grouped_lists for item in sublist], lambda a: a.predicted_gene))
    optimal = False
    num_in_list = 1 # this is the number of elements per list reurned
    best_duplicates = 0 
    splits = unique_genes_in_organism - len_unique_grouped
    #print "Begin - splits" , splits
    #print grouped_lists
    while not optimal:
        for group in itertools.combinations(grouped_lists, num_in_list):
            all_homologs_in_grouping = [item for sublist in group for item in sublist]
            #print all_homologs_in_grouping
            unique_in_set = len(MakeUnique(all_homologs_in_grouping, lambda a: a.predicted_gene))
            if unique_in_set == len_unique_grouped: # we have an optimal solution, perhaps not global optima
                duplicates =  int(math.fabs(len(all_homologs_in_grouping) - len_unique_grouped))
                if not optimal:
                    optimal = True
                    best_grouping = list(group)
                    best_duplicates = duplicates
                    best_split = splits
                elif duplicates < best_duplicates:
                    best_grouping = list(group)
                    best_duplicates = duplicates
        splits+=1
        num_in_list+=1
    #print "splits " , splits, ": best_split ", best_split
    #print "Best grouping as found by the program\n", best_grouping
    return best_grouping, best_split, best_duplicates, len_unique_grouped


def return_file_list(infolder, filter_file):
    if filter_file == '':
        return return_recursive_dir_files(infolder)   
    else:
        filter_list = [i.strip() for i in open(filter_file)]
        return [i for i in return_recursive_dir_files(infolder) if os.path.basename(i).split('.')[0] in filter_list]

# This function will take a BLAST tabular result, and remove any hits that are worse than the eval threshold provided
# The return will be a list of those hits as homolog objects.
def filter_eval(fname, e_val):
    # make a list of homolog class objects
    h_list = [Homolog.from_blast(i) for i in open(fname).readlines()]
    
    result = list(filter((lambda x: x.e_val() <= e_val), h_list))

    return result

def resolve_multiple_ORF_hits(hlist):
    result = [hlist[0]]
    
    curr_homolog = hlist[0]
    for next_homolog in hlist[1:]:
        # we have multiple BLAST hits that share a ORF, resolve using e_val
        if curr_homolog.start() == next_homolog.start():
            # If current homolog has better eval 
            if curr_homolog.e_val() <= next_homolog.e_val():
                print curr_homolog.organism(), curr_homolog.locus(), "is duplicated"
                pass
            # The current homolog has a worse eval, remove for the better example
            else:
                result.pop(-1)
                result.append(next_homolog)
                print "This totally worked", next_homolog.organism()
        else:
            result.append(next_homolog)
        # Now that we are done testing the current and next homolog against    
        curr_homolog = next_homolog

    return result


# The purpose of this function is to take a list of homologs, that have been e_val (or potenitally another means) filtered.
# The return is all homologs from organisms that contain at least one neighborhood defined by max_gap.
def return_valid_organism_homologs(hlog_list, max_gap):
    result = []
    
    org_dict = {}
    
    # Stage 1:  read the list of homologs in, and organize based on accession.  Each accession will have a list of homologs for a given operon.
    # Prior to this, the program does not sort the output.
    # This section has been tested and validated to the best of my abilities.
    #print len(hlog_list)
    for item in hlog_list:
        accession = item.accession()
        #print accession
        if accession in org_dict.keys():
            org_dict[accession].append(item)
        else:
            org_dict.update({accession:[item]})
    
    
    # Stage 2: Sort the list of homologs for each organism.  Determine gene blocks based on the max_gap criterion, 
    # and reject organisms without a gene block.  
    # This section has been tested, but not extensively.  I have added resolve_multiple_ORF_hits which is untested.
    for accession in sorted(org_dict.keys()):
        h_list = org_dict.pop(accession)
        h_list.sort(key=lambda x: x.start())
        
        # Here is where the code dealing explicitly with multiple hits to a single ORF goes:
        # currently, we only use best hit. Other resolution schemes can be envisioned.
        ORF_filetered_hlist = resolve_multiple_ORF_hits(h_list)
        
        print ORF_filetered_hlist[0].organism(), len(ORF_filetered_hlist), len(h_list)
        
        #org_dict.update({accession:h_list})
        org_dict.update({accession:ORF_filetered_hlist})
    
    # Stage 3: Organize the homologs into neighborhoods.  We will not remove any organisms at this time, we will do this in
    # the next step where appropriate.
    # This is not completely tested, but appears to be working when tested against known cases.
    
    neighborhood_dict = {}
    # This is for debug
    
    
    
    
    for accession in sorted(org_dict.keys()):
        hlist = org_dict.pop(accession)
        gene_block_list, neighborhood_found = group_homologs(hlist, max_gap)
        
        if neighborhood_found:
            neighborhood_dict.update({accession:gene_block_list})
            org_dict.update({accession:hlist})
            
        else: # do nothing, there are no neighborhoods that have been recovered
            #print "accession", accession, "is missing."
            print "Organism ", hlist[0].organism(), "is missing."
            pass
    #print sorted(neighborhood_dict.keys())
        
        
        #the following line needs to be removed, debugging purposes only
        #org_dict.update({accession:h_list})
        #if accession == 'NC_002696':
        #    print "gene_block_list", gene_block_list
    
    # Stage 4: 
    
    
    
    '''
    print org_dict.keys()
    #NC_002696
    #for h_log in org_dict['NC_000913']:
    for h_log in org_dict['NC_002696']:
        print h_log.blast_annatation(), h_log.start()
    print '\t'.join([i.blast_annatation() for i in org_dict['NC_002696']])
    
        #print '\t'.join([i.blast_annatation() for i in org_dict['NC_000913']])
    '''
    
    return result

# I think this version is more readable than those i have made in the past. 
# It can take either a sorted, or unsorted list of homologs.   
def group_homologs(lst_homologs, max_gap):
    # The gene_block_list contains a list of lists.  each list element contains genes that are no more than max_gap away.
    gene_block_list = []
    gene_block = []
    ungrouped = True
    neighborhood_found = False
    
    # bypassing pass by reference in python that leads to potentially undesirable behavior downstream
    list_homologs = [i for i in lst_homologs]
    
    # Step 1: Sort the input list by the start position of the ORF
    list_homologs.sort(key=lambda x: x.start())
    
    #for item in list_homologs:
    #    print item.accession(), item.start()
    
    # Step 2: Group homologs into gene blocks as defined my max_gap, and report these groups.
    for i in range(0, len(list_homologs)-1):
        #look at current
        start = list_homologs[i].start()
        stop = list_homologs[i].stop()
        # look at next
        start_n = list_homologs[i+1].start()
        stop_n = list_homologs[i+1].stop()
        
        if math.fabs(start - stop_n) < max_gap or math.fabs(stop - start_n) < max_gap:
            neighborhood_found = True
            if ungrouped:
                gene_block.append(list_homologs[i])
                gene_block.append(list_homologs[i+1])
                ungrouped = False
            else:
                gene_block.append(list_homologs[i+1])
        else: # get ready for next possible set of matches (in a segmented OTU)
            ungrouped = True
            if len(gene_block) > 0:
                gene_block_list.append(gene_block)
                gene_block = []
            else:
                gene_block_list.append([list_homologs[i]])
    if ungrouped:
        gene_block_list.append([list_homologs[len(list_homologs)-1]])
    else:
        gene_block_list.append(gene_block)
        
    return gene_block_list, neighborhood_found

   
def main():
    start = time.time()
    
    parsed_args = parser_code()
    
    infolder, outfolder, filter_file, num_proc, e_val, max_gap = check_options(parsed_args)
    
    print infolder, outfolder, filter_file, num_proc, e_val, max_gap
    
    file_list = return_file_list(infolder, filter_file)
    
    print file_list

    '''
    # lets test evals
    for fname in file_list:
        print "File size before", len(open(fname).readlines())
        print "File size after", len(filter_eval(fname, e_val))
    '''
    
    #'''
    for fname in file_list:
        print fname
        hlog_list = filter_eval(fname, e_val)
        return_valid_organism_homologs(hlog_list, max_gap)
    #'''
    
    # debug
    #hlog_list = filter_eval('./blast_result/NC_000913.txt', e_val)
    #return_valid_organism_homologs(hlog_list, max_gap)
    
    
    
    print time.time() - start

    # ./filter_operon_blast_results.py -f phylo_order.txt

if __name__ == '__main__':
    main()
