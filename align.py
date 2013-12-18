from utils import Bunch
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Seq import Seq

def align(seq1, seq2):
	matrix = matlist.pam60
	gap_open = -5
	gap_extend = -1
 
	alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
 
	top_aln = alns[0]
	aln_seq1, aln_seq2, score, begin, end = top_aln
	return aln_seq1, aln_seq2


def get_triplet(seq, pos, protein):
	if protein != '-':
		triplet_nucl = seq[pos : pos + 3]
		pos += 3
	else:
		triplet_nucl = '---'

	return triplet_nucl


def protein_align(seq1, seq2):
	aseq1 = aseq2 = ''
	apseq1, apseq2 = align(
		seq1.aseq[:-1],
		seq2.aseq[:-1]
	)

	mismatch = pos1 = pos2 = 0

	for n in xrange(len(apseq1)):
		triplet_nucl1 = get_triplet(seq1.seq, pos1, apseq1[n])
		pos1 += 3 if triplet_nucl1 != '---' else 0
		aseq1 += triplet_nucl1

		triplet_nucl2 = get_triplet(seq2.seq, pos2, apseq2[n])
		pos2 += 3 if triplet_nucl2 != '---' else 0
		aseq2 += triplet_nucl2

		each_mismatch = 0
		for i in xrange(3):
			each_mismatch += 1 if triplet_nucl1[i] != triplet_nucl2[i] else 0

		if apseq1[n] != apseq2[n]:
			each_mismatch *= 2

		mismatch += each_mismatch

	return Bunch(score=mismatch, aseq1=aseq1, aseq2=aseq2, aaseq1=apseq1, aaseq2=apseq2)


def pairwise_align(sequences_by_species):
	sequences = []
	for key, value in sequences_by_species.iteritems():
		sequences += value

	sequences_count = len(sequences)

	alignments = [[ None for i in xrange(sequences_count)] for j in xrange(sequences_count)]
	for i in xrange(sequences_count):
		print sequences[i].specie
		for j in xrange(sequences_count):
			alignments[i][j] = protein_align(sequences[i], sequences[j])

	return Bunch(sequences=sequences, alignments=alignments)
