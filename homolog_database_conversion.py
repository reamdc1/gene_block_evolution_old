#!/usr/bin/python

import time
import os
import sys
import argparse
import MySQLdb
from homolog4 import *

# Copyright(C) 2014 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This program's purpose is to convert a homolog data element into an entry in our database.
# it will manage insertion of single or lists of data that the user provides. 
# It is probably best to view this as an extension of the homolog class more than a seperate
# piece of coding. 

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description="Parse a homolog (or -m 8 BLAST result formatted by my software pipe) and save as a GBEER database.")
                
    parser.add_argument("-i", "--infile", dest="infile", default='/home/dave/Desktop/final_code_fork/intermediate_for_debug/unfiltered_operon/atpIBEFHAGDC.txt', metavar="FILE",
                help="A file that contains the information that you want to store in the GBEER format database.")
                
    parser.add_argument("-u", "--user", dest="user", default='root', metavar="USER",
                help="The user name for the GBEER database.")
                
    parser.add_argument("-p", "--pwd", dest="pwd", default='', metavar="PASSWORD",
                help="The password for the GBEER database.")
                
    parser.add_argument("-d", "--db", dest="db", default='gene_block', metavar="DATABASE",
                help="The name of the GBEER database in your installation.")
    
    # I do not know that we need this for the current program.  i will leave it in incase that i allow for batch infiles at a later time.
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Currently unsed, but will allow the manipulation of the number processors that you want this script to run on. The default is every CPU that the OS reports.")
                
    # I need to add the database/user/pass to this.  not sure how to best acomplish this for the project, but there is a need to these data.
                
    return parser.parse_args()


def check_options(parsed_args):
        
    if os.path.exists(parsed_args.infile):
        infile = parsed_args.infile
    else:
        print "The file %s does not exist." % parsed_args.infile
        sys.exit()
    
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)
        
    user = parsed_args.user
    pwd = parsed_args.pwd
    db = parsed_args.db
        

    return infile, num_proc, user, pwd, db

# This function handles the file traversal, and conversion into a list of homologs.
def file_to_homolog_list(infile):
    result = []
    
    handle = open(infile, 'r')
    for item in [i.strip() for i in handle.readlines()]:
        result.append(Homolog.from_blast(item))
    
    return result

# This function will convert the homolog into a GBEER database insert + data statement, and update the database.
 
 
 # I suck, so this is experimenting!
def test(usr, pwd, database):
    # db=_mysql.connect(host="localhost",user="joebob", passwd="moonpie",db="thangs")
    
    db = MySQLdb.connect("localhost", usr, pwd, database)
    cursor = db.cursor()
    cursor.execute("SELECT VERSION()")
    data = cursor.fetchone()
    print "data", data
    db.close()
    

def main():
    start = time.time()
    
    parsed_args = parser_code()
    
    infile, num_proc, user, pwd, db = check_options(parsed_args)
    
    print infile, num_proc, user, pwd, db
    
    homolog_list = file_to_homolog_list(infile)
    
    print "got here"
    #test(user, pwd, db)
    print "finished"
    
    
    print time.time() - start

if __name__ == '__main__':
    main()
