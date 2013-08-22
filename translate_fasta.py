#!/usr/bin/python
# filename: translate_fasta.py


'''
Simple script that translates a FASTA file of nucleotide sequences. Output is a single file containing
the translated sequences. 

Dependancies: biopython, python >=2.7
'''


import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser("Translates a FASTA file of nucleotide sequences.")
parser.add_argument('-in', dest='in_file', required=True, help="Name of the input FASTA file. Required")
parser.add_argument('-out', dest='out_file', required=True, help="Name of the ouput (translated) FASTA file. Required.")
args = parser.parse_args()


def main():

	# make input/output handles, clear the output file
	in_handle = open(args.in_file, 'r')
	out_handle = open(args.out_file, 'a')
	open(args.out_file, 'w').write('')

	# iterate through the input FASTA file
	for seq in SeqIO.parse(in_handle, 'fasta'):

		# grab sequence ID and translated sequence
		seq_id = seq.id
		translated_seq = str(seq.seq.translate())

		# write the sequence to the output file
		out_handle.write('>{0}\n{1}\n'.format(seq_id, translated_seq))



if __name__ == '__main__':
	main()

