#!/usr/bin/python
# filename: upload_monitor.py

###########################################################################
#
# Copyright (c) 2013 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT) 
#
###########################################################################


import os
import time
import glob
import argparse
from Bio import SeqIO
from subprocess import Popen

parser = argparse.ArgumentParser("Parses the output of IgBLAST into something suitable for import into a MySQL database")
parser.add_argument('-i', '--in', dest='input', required=True, help="The input file, to be split and processed in parallel. If a directory is given, all files in the directory will be iteratively processed.")
parser.add_argument('-o', '--out', dest='output', required=True, help="The output directory, which will contain JSON or tab-delimited output files.")
parser.add_argument('-t', '--temp', dest='temp_dir', default='', help="The directory in which temp files will be stored.  If the directory doesn't exist, it will be created.  Defaults to './temp_files'.")
args = parser.parse_args()


def file_progress(f):
	old_size = 0
	size = os.path.get_size(f)
	while old_size < size:
		old_size = size
		time.sleep(10)
		size = os.path.get_size(f)
	return f

def upload_progress(direc):
	completed = False
	while not completed:
		f = most_recent_file(direc)
		last_file = file_progress(f)
		time.sleep(30)
		if last_file == most_recent_file(direc):
			completed = True




def most_recent_file(d):
	files = filter(os.path.isfile, glob.glob(d + "/*"))
	files.sort(key=lambda x: os.path.getctime(x))
	return files[0]




