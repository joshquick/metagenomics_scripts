import numpy
from Bio import SeqIO
import sys
import pysam
from pstats import Stats
from cProfile import run

### cat sd_0001_PAO1_5k.sam | python get_alleles_from_sam.py positions.txt

def read_positions(vcfpositions):
	#all positions are 0-based (Pythonic)
	positions = set(int(pos.split('\n')[0])-1 for pos in open(vcfpositions, 'r'))
	#print 'Using %i positions' %len(positions)
	return positions

def init_array():
	#  1 2 3 4 5 ..
	#A 0 0 0 0 0 ..
	#C 0 0 0 0 0 ..
	#G 0 0 0 0 0 ..
	#T 0 0 0 0 0 ..
	genome_size = 8000000
	array = numpy.zeros((4,genome_size), dtype=int)
	return array	

def read_sam(positions, array):
	base_dict = {'A':0, 'C':1, 'G':2, 'T':3}
	#print positions
	alignment = pysam.Samfile('-', 'r')
	for line in alignment:
		#for pos in sorted(positions):
		#	if line.pos <= pos <= line.aend:
				#print line.pos, line.aend, line.alen
				#print pos-line.pos 
				#print line.seq[(pos-line.pos)-1]
				#print ord(line.qual[(pos-line.pos)-1])-33
		if line.is_unmapped: continue
		# ignore any reads with indels
		if len([f for f, v in line.cigar if f != 0]):
			continue

		read_positions = set(xrange(line.pos, line.aend))
		isec = positions.intersection(read_positions)
		if isec:
			overlap = [(pos, line.seq[pos-line.pos]) for pos in isec]
			if overlap:
				for each in overlap:
					if each[1] != 'N':
						array[(base_dict[each[1]], each[0])] += 1

	#print array.nonzero()
	return array

def write_alleles(sample_name, positions, array):
	int_dict = {0:'A', 1:'C', 2:'G', 3:'T'}
	frag = ''
	chastity_list = []
	skipped = 0
	for pos in sorted(positions):
		counts = tuple(array[:,pos])
		counts_sort = sorted(counts, reverse=True)
		if all(counts[0] > base for base in counts[1:4]):
			chastity_list.append(float(counts_sort[0]) / sum(counts_sort[:2]))
			frag += 'A'
		elif all(counts[1] > base for base in counts[0:1] + counts[2:4]):
			chastity_list.append(float(counts_sort[0]) / sum(counts_sort[:2]))
			frag += 'C'
		elif all(counts[2] > base for base in counts[0:2] + counts[3:4]):
			chastity_list.append(float(counts_sort[0]) / sum(counts_sort[:2]))
			frag += 'G'
		elif all(counts[3] > base for base in counts[0:3]):
			chastity_list.append(float(counts_sort[0]) / sum(counts_sort[:2]))
			frag += 'T'
		else:
			skipped += 1
			frag += '-'

	mean_chastity = (sum(chastity_list) / len(frag)) * 100.0
	print 'Sample: %s' %(sample_name)
	print 'Total positions: %i' %(len(positions))
	print 'Positions covered: %.3f %%' %(100 - ((skipped / len(positions)) * 100.0))
	print 'Mean chastity: %.3f %%' %(mean_chastity)
	if mean_chastity < 90:
		print 'CHASTITY WARNING: Mixed samples can severely affect accuracy of placement'

	id = sample_name
	file_out = open('alleles/%s.fasta' %(sample_name), 'w')
	print >>file_out, '>%s_new\n%s\n' %(id, frag),
	file_out.close()

if __name__ == '__main__':
	positions = read_positions(sys.argv[2])
	positions_set = set(positions)
	sample_name = sys.argv[1]
	array = init_array()
	run('read_sam(positions_set, array)', 'read_sam.stats')
	stats = Stats('read_sam.stats')
	write_alleles(sample_name, positions, array)
