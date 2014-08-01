#!/usr/bin/python
# filename: barcode_generator.py


import argparse
import random
import itertools

from Bio import pairwise2


parser = argparse.ArgumentParser("Generates a specified number of barcode sequences suitable for indexing on Illumina runs.")
parser.add_argument('-o', '--out', dest='output_file', default="", help="Output file name. Required.")
parser.add_argument('-b', '--prefix', dest='prefix', default="", help="Prefix for barcode names. 'ABC' would result in primers named ABC001, ABC002, etc. Default is ''.")
parser.add_argument('-n', '--num', dest='num', type=int, default=96, help="Number of barcodes to generate. Default is 96.")
parser.add_argument('-l', '--length', dest='length', type=int, default=10, help="Barcode length, in nucleotides. Default is 10.")
parser.add_argument('-m', '--mismatch', dest='mismatch', type=int, default=3, help="Minimum number of mismatches between barcodes. Default is 3.")
parser.add_argument('-s', '--seed', dest='seed', type=str, default='1234', help="Seed for the random barcode generator. Default is '1234'.")
parser.add_argument('-x', '--nextseq', dest='nextseq', default=False, action='store_true', help="Makes NextSeq-compatible barcodes (at least one non-G residue in first five positions).")
args = parser.parse_args()

def all_possibilities(length):
	nucs = ['A', 'C', 'G', 'T']
	return [''.join(p) for p in itertools.product(nucs, repeat=length)]

def distance(b, barcodes):
	scores = []
	for barcode in barcodes:
		scores.append(pairwise2.align.globalxx(b, barcode, one_alignment_only=True, score_only=True))
	if max(scores) <= args.length - args.mismatch:
		return True
	return False

def filter_mismatches(barcodes):
	random.shuffle(barcodes)
	filtered = [barcodes[0],]
	for b in barcodes[1:]:
		if distance(b, filtered):
			filtered.append(b)
			if len(filtered) % 25 == 0:
				print len(filtered)
		if len(filtered) == args.num:
			break
	print len(filtered)
	return filtered

def filter_nextseq(barcodes):
	return [b for b in barcodes if b[:5] != 'GGGGG']

def calc_stats(barcodes):
	print '\nNucleotide Frequencies:'
	print 'Pos\tA\tC\tG\tT'
	residues = {}
	for i in range(args.length):
		residues[i] = {}
		for b in barcodes:
			residues[i][b[i]] = residues[i][b[i]] + 1 if b[i] in residues[i] else 1
	for i in range(args.length):
		print '{}\t'.format(i+1),
		for n in ['A', 'C', 'G', 'T']:
			count = residues[i][n] if n in residues[i] else 0
			print '{}\t'.format(count),
		print '' 

def write_output(barcodes):
	chars = len(str(args.num))
	fasta_string = ''
	for i, b in enumerate(barcodes):
		zeroes = chars - len(str(i+1))
		fasta_string += '{}{}{}\t{}\n'.format(args.prefix, '0'*zeroes, str(i+1), b)
	open(args.output_file, 'w').write(fasta_string)


def main():
	barcodes = all_possibilities(args.length)
	filtered_barcodes = filter_mismatches(barcodes)
	if args.nextseq:
		filtered_barcodes = filter_nextseq(filtered_barcodes)
	calc_stats(filtered_barcodes)
	write_output(filtered_barcodes)



if __name__ == '__main__':
	random.seed(args.seed)
	main()


