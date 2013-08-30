#!/usr/local/bin/python
# filename: fastq_to_fasta.py

import argparse

parser = argparse.ArgumentParser("Converts a file of sequences in FASTQ format to FASTA file\n")
parser.add_argument('-in', dest='seq_file',required=True, help="The input file, in FASTQ format")
parser.add_argument('-out', dest='output_file', required=True, help="The output file, in FASTA format. Defaults to \"<input>.fasta\"")
args = parser.parse_args()


def main():

	# clear the output file
	open(args.output_file, 'w').write('')
	 
	input_handle = open(args.seq_file, "rU")
	output_handle = open(args.output_file, "a")


	# iterate through the file and grab the sequence IDs and sequences
	next_is_sequence = False
	for line in input_handle.read().split('\n'):

		print line

		# for blank lines
		if len(line) < 1:
			continue

		# am I expecting a sequence on this line?
		if next_is_sequence:
			output_handle.write(line + '\n')
			next_is_sequence = False
			continue

		# if the line is a sequence ID line
		if line[0] == '@':
			output_handle.write('>' + line[1:] + '\n')
			next_is_sequence = True
			continue
	 
	output_handle.close()
	input_handle.close()


if __name__ == '__main__':
	main()