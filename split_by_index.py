#! /usr/bin/python
# filename: split_by_index.py


import argparse
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


parser = argparse.ArgumentParser("Splits a file of sequences based on indexes at or near the start of the read.")
parser.add_argument('-in', dest='input', required=True, help="The input file, contained indexed sequences in either FASTA or FASTQ format.")
parser.add_argument('-out', dest='out_prefix', default="", help="Prefix for the output files.	Defaults to '<input>_<index_seq>.fasta'.")
parser.add_argument('-indexes', dest='indexes', required=True, help="A fasta file containing the index sequences. All indexes must be the same length.")
parser.add_argument('-format', dest='input_format', default='', help="Define the input file format (either fasta or fastq).  Default is fasta.")
parser.add_argument('-offset', dest='offset', default=0, type=int, help="Number of bases to ignore at the start of the sequence when searching for the index. Useful if a defined number of ambiguous nucleotides were inserted at the start of the sequencing read (to aid cluster identification).")
args = parser.parse_args()




def clear_output_files(indexes, out_prefix):

	# clear each of the demultiplexed output files
	for index in indexes:
		indexed_file = '{0}_{1}.fasta'.format(out_prefix, index)
		open(indexed_file, 'w').write('')

	# clear the 'unassigned' output file
	open(out_prefix + '_unassigned.fasta', 'w').write('')



def parse_indexes():

	# set up some variables
	index_handle = open(args.indexes, 'r')
	indexes = []

	# parse the index file
	for record in SeqIO.parse(index_handle, 'fasta'):
		indexes.append(str(record.seq))

	# get the index length
	index_length = len(indexes[0])

	return indexes, index_length



def get_index(seq, index_length):

	return seq[args.offset : args.offset + index_length]



def main():

	# define the input handle
	input_handle = open(args.input, 'r')

	# parse the index file
	indexes, index_length = parse_indexes()

	# clear the output files
	clear_output_files(indexes, args.out_prefix)

	# build a dict to hold the output handles
	out_handles = {}
	for index in indexes:
		out_file = '{0}_{1}.fasta'.format(args.out_prefix, index)
		out_handles[index] = open(out_file, 'a')

	# add an 'unassigned' handle to the out_handles dict
	out_handles['unassigned'] = open(args.out_prefix + '_unassigned.fasta', 'a')

	# iterate through each sequence
	for record in SeqIO.parse(input_handle, args.format):

		# identify the index
		seq_id = record.id
		sequence = str(record.seq)
		index = get_index(sequence, index_length)

		# determine whether the found index is an appropriate index
		if index in indexes:

			# write the sequence to the appropriate output file
			out_handles[index].write('>{0}\n{1}\n'.format(seq_id, sequence[args.offset + index_length:]))

		# if the index isn't in the indexes list
		else:

			# write the sequence to the 'unassigned' output file
			out_handles['unassigned'].write('>{0}\n{1}\n'.format(seq_id, sequence[args.offset + index_length:]))


if __name__ == '__main__':
	main()

