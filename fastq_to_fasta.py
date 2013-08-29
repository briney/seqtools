#!/usr/local/bin/python
# filename: fastq_to_fasta.py


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

parser = argparse.ArgumentParser("Converts a file of sequences in FASTQ format to FASTA file\n")
parser.add_argument('-in', dest='seq_file',required=True, help="The input file, in FASTQ format")
parser.add_argument('-out', dest='output_file', default="",help="The output file, in FASTA format. Defaults to \"<input>.fasta\"")
parser.add_argument('-sra', dest='from_sra', action='store_true', default=False, help="Use this flag if parsing sequences from the SRA archive, which have screwed up sequence IDs.")
parser.add_argument('-v', dest='verbose', action="store_true", help="Print values to screen")
args = parser.parse_args()


if args.output_file == "":
	args.output_file = args.seq_file.split(".")[0]+".fasta"
 
input_handle = open(args.seq_file, "rU")
output_handle = open(args.output_file, "w")

sequences = (rec for rec in SeqIO.parse(input_handle, "fastq"))
for seq in sequences:
	if args.from_sra:
		output_handle.write('>%s\n%s\n' %(seq.description.split()[1], str(seq.seq)))
	else:
		output_handle.write('>%s\n%s\n' %(seq.id, str(seq.seq)))
 
output_handle.close()
input_handle.close()