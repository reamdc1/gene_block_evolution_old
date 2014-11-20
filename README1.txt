This document is intended to grant users of this software a fighting chance for susccessful use.
The general organization will be the in order documentation of use for the entire program pipeline.

First things first.  The entire software package that this document relates to is released under GPL version 3.
The full licence can be located here: http://www.gnu.org/licenses/lgpl.html


This can be run at any time, but really does not belong in any formal pipeline for the main program, so it comes 
First.
0) regulondb_dl_parse.py

Download and parse a regulonDB operon file, then determine information about
the component genes by investigating the refrence organisms.

optional arguments:
  -h, --help            show this help message and exit
  -i FOLDER, --infolder FOLDER
                        Folder containing all genbank files for use by the
                        program.
  -o FOLDER, --outfolder FOLDER
                        Folder where results will be stored.
  -f FILE, --filter FILE
                        File restrictiong which accession numbers this script
                        will process. If no file is provided, filtering is not
                        performed.
  -n INT, --num_proc INT
                        Number of processors that you want this script to run
                        on. The default is every CPU that the system has.
  -d, --download        Add this option if you wish to download the regulonDB
                        operon, otherwise the program will assume that you
                        have already done this step.
  -u URL, --url URL     A file that contains the BLAST query for every gene of
                        interest in the dataset.
  -e, --experimantal    Add this option if you wish to download the regulonDB
                        operon, otherwise the program will assume that you
                        have already done this step.
  -m INT, --min_genes INT
                        Minum number of genes that an operon must contain
                        before it can be considered for further analysis. The
                        default is 5 because that is what we are currently
                        using in the study.


Example:
./regulondb_dl_parse.py

1) file_lists.py

Determine information about the genbank files under study, and report this
information for use by other modules.

optional arguments:
  -h, --help            show this help message and exit
  -i FOLDER, --infolder FOLDER
                        Folder containing all genbank files for use by the
                        program.
  -o FOLDER, --outfolder FOLDER
                        Folder where results will be stored.
  -f FILE, --filter FILE
                        File restrictiong which accession numbers this script
                        will process. If no file is provided, filtering is not
                        performed.
  -n INT, --num_proc INT
                        Number of processors that you want this script to run
                        on. The default is every CPU that the system has.

Examples:
./file_lists.py -i /home/dave/Desktop/all_genbank -o genbank_pathway_lists_proteobacteria -f phylo_order.txt
./file_lists.py -o genbank_pathway_lists_proteobacteria -f phylo_order.txt


2) format_db.py

Convert all genbank files found in a specified folder, and optionally a file
containing the accessions that you wish to include, and create BLAST
searchable databases from them.

optional arguments:
  -h, --help            show this help message and exit
  -i FOLDER, --infolder FOLDER
                        Folder containing all genbank files for use by the
                        program.
  -o FOLDER, --outfolder FOLDER
                        Folder where the BLAST searchable databases will be
                        stored.
  -f FILE, --filter FILE
                        File restrictiong which accession numbers this script
                        will process. If no file is provided, filtering is not
                        performed.
  -n INT, --num_proc INT
                        Number of processors that you want this script to run
                        on. The default is every CPU that the system has.
  -p, --protein         Flag to toggle mode from protein to DNA sequences for
                        use as database construction. The default is protein.


Examples:
./format_db.py -f ./phylo_order.txt
./format_db.py -i /home/dave/Desktop/all_genbank -o ./db1/ -f ./phylo_order.txt



3) 



