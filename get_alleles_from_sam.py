import numpy
from Bio import SeqIO
import sys
import pysam
from collections import OrderedDict
#import vcf
#from cProfile import run
#from pstats import Stats

USE_CHASTITY = True

### cat sam_file | python get_alleles_from_sam.py sample_name positions.txt file_out

def read_positions(positions_txt):
	# all positions are 0-based (Pythonic)
	positions = OrderedDict()
	for line in open(positions_txt, 'r'):
		cols = line.strip().split(',')
		if cols[0] not in positions.keys():
			positions[cols[0]] = set()
			positions[cols[0]].add(int(cols[1]) - 1)
		else:
			positions[cols[0]].add(int(cols[1]) - 1)
	return positions

def init_array(refs):
	#  1 2 3 4 5 ...
	#A 0 0 0 0 0 ...
	#C 0 0 0 0 0 ...
	#G 0 0 0 0 0 ...
	#T 0 0 0 0 0 ...
	array = {}
	positions = read_positions(sys.argv[2])
	for ref in positions.keys():
		array[ref] = numpy.zeros((4, max(positions[ref]) + 1), dtype=int)
	return array, positions

def read_sam():
	base_dict = {'A':0, 'C':1, 'G':2, 'T':3}
	alignment = pysam.Samfile('-', 'r')
	array, positions = init_array(alignment.references)
	for line in alignment:
		# ignore any unmapped reads
		if line.is_unmapped:
			continue
		# ignore any reads with indels
		#if len([op for op, length in line.cigar if op != 0]):
			#check for isec first then check length of indel?
		chrom = alignment.getrname(line.tid)
		read_positions = set(xrange(line.pos, line.aend))
		try:
			isecs = positions[chrom].intersection(read_positions)
		except KeyError:
			continue
		if isecs:
			#overlap = [(pos, line.seq[pos-line.pos]) for pos in isec]
			#quality = [(pos, ord(line.qual[pos-line.pos])-33) for pos in isec]
			aligned_pairs = dict((ref, query) for (query, ref) in line.get_aligned_pairs())
			for isec in isecs:
				if aligned_pairs[isec]:
					read_base = line.seq[aligned_pairs[isec]]
					if read_base != 'N':
							array[chrom][(base_dict[read_base], isec)] += 1
	return array, positions

def write_alleles():
	sample_name = sys.argv[1]
	array, positions = read_sam()
	frag = ''
	chastity_list = []
	gaps = 0
	for chrom in positions.keys():
		for pos in sorted(positions[chrom]):
			counts = tuple(array[chrom][:,pos])
			counts_sort = sorted(counts, reverse=True)
			if all(counts[0] > base for base in counts[1:4]):
				# chastity is greatest / (greatest + second greatest)
				if USE_CHASTITY:
					chastity_list.append(float(counts_sort[0]) / 
						     sum(counts_sort[:2]))
				frag += 'A'
			elif all(counts[1] > base for base in counts[0:1] + counts[2:4]):
				if USE_CHASTITY:
					chastity_list.append(float(counts_sort[0]) / 
						     sum(counts_sort[:2]))
				frag += 'C'
			elif all(counts[2] > base for base in counts[0:2] + counts[3:4]):
				if USE_CHASTITY:
					chastity_list.append(float(counts_sort[0]) / 
						     sum(counts_sort[:2]))
				frag += 'G'
			elif all(counts[3] > base for base in counts[0:3]):
				if USE_CHASTITY:
					chastity_list.append(float(counts_sort[0]) / 
						     sum(counts_sort[:2]))
				frag += 'T'
			else:
				gaps += 1
				frag += '-'

	total = sum([len(positions[chrom]) for chrom in positions.keys()])
	print 'Sample: %s' %(sample_name)
	print 'Total positions: %i' %(total)
	print 'Gaps: %i' %(gaps)
	print 'Positions covered: %.2f %%' %(100 - (float(gaps) / total * 100.0))
	print 'Sample: %s' %(sample_name)
	if USE_CHASTITY:
		mean_chastity = (sum(chastity_list) / (len(frag) - gaps)) * 100.0
		print 'Mean chastity: %.2f %%' %(mean_chastity)
		if mean_chastity < 90:
			print 'CHASTITY WARNING: Mixed samples can severely affect accuracy of placement'

	id = sample_name
	with open(sys.argv[4], 'w') as file_out:
		print >>file_out, '>%s_new\n%s\n' %(id, frag),

if __name__ == '__main__':
	#run('write_alleles()', 'stats')
	#stats = Stats('stats')
	write_alleles()
