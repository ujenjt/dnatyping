import re
import sys
from utils import Bunch

def sequences(stream):
	sequence = ''
	not_eof = True
	while not_eof:
		pos = stream.tell()
		newline = stream.readline()

		if newline:
			description = newline if newline[0] == '>' and not sequence else description

		if newline == '':
			not_eof = False
			yield (description, sequence)
		elif newline[0] == '>' and sequence:
			stream.seek(pos)
			yield (description, sequence)
			sequence = ''
		else:
			sequence += newline if newline[0] != '>' else ''


def translate_dna(sequence):
    gencode = {
	    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
	    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    aseq = ''
    for n in range(0,len(sequence),3):
        if sequence[n:n+3] in gencode:
            aseq += gencode[sequence[n:n+3]]

    return aseq


def parse_data(inp):
	merged_sequences = dict()

	for s in sequences(inp):
		seq = s[1].replace('\n', '').replace('\r', '')
		matchObj = re.match(r'.*(Vibrio [^\s.]*).*', s[0], re.M|re.I)
		specie = matchObj.group(1)

		if specie != 'Vibrio sp' and specie != 'Vibrio genomosp':
			sequence = merged_sequences.get(seq, Bunch(seq=seq, specie=specie, merged_count=0, row_headers=[], aseq=translate_dna(seq)))
			sequence.merged_count += 1
			sequence.row_headers.append(s[0])

			merged_sequences[seq] = sequence

	return merged_sequences


def import_data(inp):
	merged_sequences = parse_data(inp)

	sequences_by_species = dict()
	for key, value in merged_sequences.items():
		sequences = sequences_by_species.get(value.specie, [])
		sequences.append(value)
		sequences_by_species[value.specie] = sequences

	return sequences_by_species
