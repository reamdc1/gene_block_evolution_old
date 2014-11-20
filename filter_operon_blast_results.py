#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import sys
import argparse
from homolog4 import *
from collections import defaultdict
import itertools

# Copyright(C) 2014 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description="This program will be used to remove spurious results from a BLAST search organized by operon.  The result will be the basis of the database for GBEER and for other anaylsis.")
                
    parser.add_argument("-i", "--infolder", dest="infolder", default='./blast_result/', metavar="FOLDER",
                help="A folder that contains the operon BLAST results.")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="FOLDER", default='./optimized_operon/',
                help="Folder where the filtered results will be stored. Default is the folder './optimized_operon/'.")

    parser.add_argument("-f", "--filter", dest="filter", default='', metavar="FILE",
                help="A file that contains the accession numbers of the organisms that are under investigation.  All others will be omitted from analysis an results.")            
    
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
                
    parser.add_argument("-e", "--eval", dest="eval", default='1e-10', metavar="FLOAT", type=float,
                help="Use this option to change the eval for the BLAST search that is permitted. Useful if you would like to investigate what altering the eval threshold will do to your results.")
                
    parser.add_argument("-g", "--max_gap", dest="maw_gap", metavar="INT", default = 500, type=int,
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

# this function will take a list of homologs and return a list of lits
# each element will be a list of grouped homologs. the criteria is that
# each item will be no further than INTERGENIC_MAX_LEN away from its closest neighbor
# strand is completely ignored for this function, though i do not know that is best.
def GroupHomologs(list_homologs):
    comprehensive_list = []
    local_group = []
    ungrouped = True
    for i in range(0, len(list_homologs)-1):
        #look at current
        start = list_homologs[i].start
        stop = list_homologs[i].stop
        # look at next
        start_n = list_homologs[i+1].start
        stop_n = list_homologs[i+1].stop
        
        if math.fabs(start - stop_n) < INTERGENIC_MAX_LEN or math.fabs(stop - start_n) < INTERGENIC_MAX_LEN:
            if ungrouped:
                local_group.append(list_homologs[i])
                local_group.append(list_homologs[i+1])
                ungrouped = False
            else:
                local_group.append(list_homologs[i+1])
        else: # get ready for next possible set of matches (in a segmented OTU)
            ungrouped = True
            if len(local_group) > 0:
                comprehensive_list.append(local_group)
                local_group = []
            else:
                comprehensive_list.append([list_homologs[i]])
    if ungrouped:
        comprehensive_list.append([list_homologs[len(list_homologs)-1]])
    else:
        comprehensive_list.append(local_group)
    return comprehensive_list, [[j for j in i] for i in comprehensive_list if len(i) > 1] # this last one is only grouped homologs

def LongestGroup(lst):
    if len(lst) == 0:
        return 0
    else:
        return max([len(i) for i in lst])



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


# The purpose of this function is to take a list of homologs, that have been e_val (or potenitally another means) filtered.
# The return is all homologs from organisms that contain at least one neighborhood defined by max_gap.
def return_valid_organism_homologs():
    result = []
    
    org_dict = {}
    
    
    return result
    
    


   
def main():
    start = time.time()
    
    parsed_args = parser_code()
    
    infolder, outfolder, filter_file, num_proc, e_val, max_gap = check_options(parsed_args)
    
    print infolder, outfolder, filter_file, num_proc, e_val, max_gap
    
    file_list = return_file_list(infolder, filter_file)
    #print len(file_list)
    
    # lets test evals
    
    for fname in file_list:
        print "File size before", len(open(fname).readlines())
        print "File size after", len(filter_eval(fname, e_val))
    
    print time.time() - start

    # ./filter_operon_blast_results.py -f phylo_order.txt

if __name__ == '__main__':
    main()
