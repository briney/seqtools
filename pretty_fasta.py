#!/usr/bin/python
# filename: pretty_fasta.py


'''
Simple script that converts an 'ugly' FASTA file -- those for which the sequence contains line 
breaks -- into a 'pretty' FASTA file with each sequence on a single line.  

Dependancies: biopython, python >=2.7
'''


import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser("Converts an 'ugly' FASTA file (with sequences spanning multiple lines) to a 'pretty' one.")
parser.add_argument('-in', dest='in_file', required=True, help="Name of the input FASTA file. Required")
parser.add_argument('-out', dest='out_file', required=True, help="Name of the ouput FASTA file. Required.")
args = parser.parse_args()


def main():

	# make input/output handles, clear the output file
	in_handle = open(args.in_file, 'r')
	out_handle = open(args.out_file, 'a')
	open(args.out_file, 'w').write('')

	# iterate through the input FASTA file
	for seq in SeqIO.parse(in_handle, 'fasta'):

		# grab the sequence and the sequence ID
		seq_seq = str(seq.seq)
		seq_id = seq.id

		# write the sequence to the output file
		out_handle.write('>{0}\n{1}\n'.format(seq_id, seq_seq))



if __name__ == '__main__':
	main()

