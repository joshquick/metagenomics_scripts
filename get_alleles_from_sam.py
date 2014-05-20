import numpy
from Bio import SeqIO
import sys
import pysam
from pstats import Stats
from cProfile import run

### cat sd_0001_PAO1_5k.sam | python get_alleles_from_sam.py sample_name positions.txt

def read_positions(positions_file):
	# all positions are 0-based (Pythonic)
	positions = set(int(pos.split('\n')[0])-1 for pos in open(positions_file, 'r'))
	return positions

def init_array(positions):
	#  1 2 3 4 5 ..
	#A 0 0 0 0 0 ..
	#C 0 0 0 0 0 ..
	#G 0 0 0 0 0 ..
	#T 0 0 0 0 0 ..
	genome_size = max(positions) + 1
	array = numpy.zeros((4,genome_size), dtype=int)
	return array

def read_sam(positions, array):
	base_dict = {'A':0, 'C':1, 'G':2, 'T':3}
	alignment = pysam.Samfile('-', 'r')
	for line in alignment:
		# ignore any unmapped reads
		if line.is_unmapped: continue
		# ignore any reads with indels
		if len([f for f, v in line.cigar if f != 0]):
			continue

		read_positions = set(xrange(line.pos, line.aend))
		isec = positions.intersection(read_positions)
		if isec:
			overlap = [(pos, line.seq[pos-line.pos]) for pos in isec]
			#quality = [(pos, ord(line.qual[pos-line.pos])-33) for pos in isec]
			if overlap:
				for each in overlap:
					if each[1] != 'N':
						array[(base_dict[each[1]], each[0])] += 1
	return array

def write_alleles(sample_name, positions, array):
	frag = ''
	chastity_list = []
	gaps = 0
	for pos in sorted(positions):
		counts = tuple(array[:,pos])
		counts_sort = sorted(counts, reverse=True)
		if all(counts[0] > base for base in counts[1:4]):
			# chastity is greatest / (greatest + second greatest)
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
			gaps += 1
			frag += '-'

	mean_chastity = (sum(chastity_list) / len(frag)) * 100.0
	print 'Sample: %s' %(sample_name)
	print 'Total positions: %i' %(len(positions))
	print 'Gaps: %i' %(gaps)
	print 'Skipped: %i' %(skipped)
	print 'Positions covered: %.3f %%' %(100 - (float(skipped) / len(positions) * 100.0))
	print 'Mean chastity: %.3f %%' %(mean_chastity)
	if mean_chastity < 90:
		print 'CHASTITY WARNING: Mixed samples can severely affect accuracy of placement'

	id = sample_name
	with open('alleles/%s.fasta' %(sample_name), 'w') as file_out:
		print >>file_out, '>%s_new\n%s\n' %(id, frag),

if __name__ == '__main__':
	positions = read_positions(sys.argv[2])
	sample_name = sys.argv[1]
	array = init_array(positions)
	read_sam(positions, array)
	write_alleles(sample_name, positions, array)
