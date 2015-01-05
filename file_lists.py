#!/usr/bin/python

# Copyright(C) 2014 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

import os
import argparse
import time
import sys
import re
import Bio
from Bio import SeqIO,SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
import simplejson as json
from multiprocessing import Pool
import cPickle as pickle
import gc


# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description='Determine information about the genbank files under study, and report this information for use by other modules.')

    parser.add_argument("-i", "--infolder", dest="infolder", metavar="FOLDER", default='/home/dave/Desktop/all_genbank',
                help="Folder containing all genbank files for use by the program.")
                 
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="FOLDER", default='./genbank_pathway_lists/',
                help="Folder where results will be stored.")
    
    parser.add_argument("-f", "--filter", dest="filter", metavar="FILE", default='',
                help="File restrictiong which accession numbers this script will process. If no file is provided, filtering is not performed.")
                
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
                
    return parser.parse_args()
    
    
# A function that will check the command line input for errors. If serious errors exist, it will exit.
def check_options(parsed_args):

    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder    
    else:
        print "The folder %s does not exist." % parsed_args.infolder
        sys.exit()
        
    if infolder[-1] != '/':
        infolder = infolder + '/'
    
    # if the directory that the user specifies does not exist, then the program makes it for them. 
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    if parsed_args.outfolder[-1] != '/':
       outfolder = parsed_args.outfolder + '/'
    else:
       outfolder = parsed_args.outfolder 
    
    if parsed_args.filter == '' or os.path.exists(parsed_args.filter):
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

    return infolder, outfolder, filter_file, num_proc


#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# This function will take in the folder where the genomes reside.  By deafult i have this under Desktop/all_genbank/
# but this could be different. This function will make a master file list, genome only file list, master organism file list,
# and a filtered set. Filtering will be at the NC number level right now. I am not allowing renaming of the file names, because
# there is really no reason for this nonsense.
def return_genbank_paths_of_interest(in_folder, out_folder, filter_fname):
    
    do_filter = not filter_fname == ''
    
    if do_filter:
        filter_list = [i.strip() for i in open(filter_fname).readlines()] 
    else:
        filter_list = [] 
    flist = returnRecursiveDirFiles(in_folder)
    
    # This first block is just the raw dump of a recursive folder walk. 
    # This code is totally irrelevant, but i will keep anyway, at least for now.
    handle = open(out_folder+"all_genbank_paths.txt", 'w')
    handle.write('\n'.join(flist))
    handle.close()
    
    # This block of code is to provide a list of chromosomes and plasmids that are contained within an organisms' folder
    # The formatting will be organism name (as read from the folder name, and generally/always? corresponds to the common
    # name of the organism, then a tab delineated list of file paths.
    f_dict = {}
    for fname in flist:
        folder = fname.split('/')[-2]
        if folder not in f_dict.keys():
            f_dict.update({folder:[fname]})
        else:
            f_dict[folder].append(fname)

    handle = open(out_folder+"grouped_by_organism_all_genbank_paths.txt", 'w')
    for org in sorted(f_dict.keys()):
        lst = [org] + sorted(f_dict[org])
        handle.write('\t'.join(lst)+'\n')
    handle.close()

    res = []
    if do_filter:
        for fname in flist:
            if fname.split('/')[-1].split('.')[0] in filter_list:
                res.append(fname)
    else:
        for fname in flist:
            res.append(fname)
            
    handle = open(out_folder+"filtered_list_genbank_paths.txt", 'w')
    handle.write('\n'.join(res))
    handle.close()
    
    return res
    
def convert_genbank(genbank_path):

    handle = open(genbank_path)
    seq_record = SeqIO.parse(handle, "genbank").next()
    

    nc_number = seq_record.name
    
    # genbank information about the sequence found in the file. 
    accession = seq_record.id.split('.')[0]
    comment = seq_record.annotations['comment']
    data_file_division = seq_record.annotations['data_file_division'] # this is not very useful, since it is historical and does not reflect current taxonomy
    date = seq_record.annotations['date'] # this should be the date that the file was updated last. uesful if I wish to autoupdate for more current versions of a genome
    gi = seq_record.annotations['gi'] # I do not plan on using this, but it may be handly to have, so I will include anyway
    key_words = seq_record.annotations['keywords']
    organism = seq_record.annotations['organism'].replace(' ', '_')
    sequence_version = seq_record.annotations['sequence_version'] # this is an integer
    taxonomy = seq_record.annotations['taxonomy'] # this will be a list
    source_folder = genbank_path.split('/')[len(genbank_path.split('/'))-2]

    #####################################################################################################################################
    # Testing of this code block indicates that there is some type of memory error in SeqIO.parse, so if this is used on a large number #
    # of sequences, you need to have a large memory system. (Currently the minumum for all sequences is 32 GB [Aug, 2013])              #
    #####################################################################################################################################
    
    # this block of code determines the type of sequence [chromosome, plasmid] and if there are more than one chromosome in the organism.
    # In the case of mutiple chromosomes, the code will assign a name to each, based on the information provided in the genbank file.
    description = seq_record.description
    if re.search("complete genome", description) is not None: # This is a complete genome
        seq_type = 'complete_genome'
        seq_name = 'complete_genome'
    elif re.search("hromosome", description) is not None: # we have a segmented chromosome (more than one chromosome)
        seq_type = 'chromosome'
        seq_name = 'chromosome' + description.split('hromosome')[1].split(',')[0].strip()
    elif re.search("lasmid", description) is not None: # we have a plasmid
        seq_type = 'plasmid'
        seq_name = description.split('lasmid')[1].split(',')[0].strip()
    else:
        print 'error', accession
        print description
        seq_type = 'complete_genome'
        seq_name = 'complete_genome'
        
    # This is intended to stop memory leak indicated above, will see    
    handle.close()
    gc.collect()
        
    return accession, {'seq_type': seq_type, 'seq_name': seq_name, 'taxonomy': taxonomy, 'organism': organism, 'source_folder': source_folder, 'file_path': genbank_path}

def convert_genbank2(genbank_path):

    
    #handle = open(genbank_path)
    with open(genbank_path, 'r') as handle:
        seq_record = SeqIO.parse(handle, "genbank").next()
        #handle.close()
    
        nc_number = seq_record.name
        
        # genbank information about the sequence found in the file. 
        accession = seq_record.id.split('.')[0]
        comment = seq_record.annotations['comment']
        data_file_division = seq_record.annotations['data_file_division'] # this is not very useful, since it is historical and does not reflect current taxonomy
        date = seq_record.annotations['date'] # this should be the date that the file was updated last. uesful if I wish to autoupdate for more current versions of a genome
        gi = seq_record.annotations['gi'] # I do not plan on using this, but it may be handly to have, so I will include anyway
        key_words = seq_record.annotations['keywords']
        organism = seq_record.annotations['organism'].replace(' ', '_')
        sequence_version = seq_record.annotations['sequence_version'] # this is an integer
        taxonomy = seq_record.annotations['taxonomy'] # this will be a list
        source_folder = genbank_path.split('/')[len(genbank_path.split('/'))-2]
    
        #####################################################################################################################################
        # Testing of this code block indicates that there is some type of memory error in SeqIO.parse, so if this is used on a large number #
        # of sequences, you need to have a large memory system. (Currently the minumum for all sequences is 32 GB [Aug, 2013])              #
        #####################################################################################################################################
        
        # this block of code determines the type of sequence [chromosome, plasmid] and if there are more than one chromosome in the organism.
        # In the case of mutiple chromosomes, the code will assign a name to each, based on the information provided in the genbank file.
        description = seq_record.description
        if re.search("complete genome", description) is not None: # This is a complete genome
            seq_type = 'complete_genome'
            seq_name = 'complete_genome'
        elif re.search("hromosome", description) is not None: # we have a segmented chromosome (more than one chromosome)
            seq_type = 'chromosome'
            seq_name = 'chromosome' + description.split('hromosome')[1].split(',')[0].strip()
        elif re.search("lasmid", description) is not None: # we have a plasmid
            seq_type = 'plasmid'
            seq_name = description.split('lasmid')[1].split(',')[0].strip()
        else:
            print 'error', accession
            print description
            seq_type = 'complete_genome'
            seq_name = 'complete_genome'
    
    # This is intended to stop memory leak indicated above, will see    
    handle.close()
    gc.collect()
        
    return accession, {'seq_type': seq_type, 'seq_name': seq_name, 'taxonomy': taxonomy, 'organism': organism, 'source_folder': source_folder, 'file_path': genbank_path}


def parallel_convert_genbank(genbank_path_list, num_proc, outfile = './genbank_pathway_lists/nc_information_dict.json'):

    pool = Pool(processes = num_proc)
    #result = dict(pool.map(convert_genbank, genbank_path_list))
    result = dict(pool.map(convert_genbank2, genbank_path_list))

    # Store the whole data structure about each of the genbank sequences as JSON because we will be using this for web services later.
    handle = open(outfile, 'w')
    json.dump(result, handle)
    handle.close()
    
    file_name, file_extension = os.path.splitext(outfile)
    
    pickle.dump(result, open(outfile.replace(file_extension, '.p'), "wb"))


def main():

    # Timer, used during debug to determine the fastest implementation for a code block
    start = time.time()

    parsed_args = parser_code()
    
    infolder, outfolder, filter_file, num_proc = check_options(parsed_args)

    genbank_path_list = return_genbank_paths_of_interest(infolder, outfolder, filter_file)
    
    # The program returns a dictionary containing information about every genbank file that exists in the genbank directory.
    # This operation takes up about 29GB+ memory (in Aug 2013) and runs in about 6-7 min. To have it report only about organisms
    # that are currently being investigated, run the second version of this function. It will take a lot less time/memory.
    # however, there is no detrimental consequence ot keeping things how they are right now.
    
    # New code runes it in under 25GB and 6:20 (Apr 2014). Memory deallocation is still a problem though.
    
    # Newer code (Oct, 2014) can run on 16GB quad core in 11 min. Memory deallocation is still a problem
    
    parallel_convert_genbank(genbank_path_list, num_proc)

    print time.time() - start
    
    # sample cmds that will run using default settings
    #./file_lists.py
    # something a bit better
    # ./file_lists.py -o genbank_pathway_lists_proteobacteria -f phylo_order.txt
    #or
    # ./file_lists.py -i /home/dave/Desktop/all_genbank -o genbank_pathway_lists_proteobacteria -f phylo_order.txt
    

if __name__ == '__main__':
    main()
