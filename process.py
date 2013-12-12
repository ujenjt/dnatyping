from parse import import_data
from align import pairwise_align
from present import present

inp = open('input_data.txt', 'r')

sequences_by_species = import_data(inp)
alignment_data = pairwise_align(sequences_by_species)
present(alignment_data)