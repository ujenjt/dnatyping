from string import Template
from utils import Bunch
import pystache
import shutil
import os.path

tmpl_file = open('assets/template.html', 'r')
output_file = open('out/result.html', 'w')

def generate_template(data, names):
	tmpl = Template(tmpl_file.read())
	vals = { 
		'table': template_table_data(append_names(convert(names), data)), 
		'head': template_table_data([ append_x(convert(names)) ])
	}
	output_file.write(tmpl.substitute(vals))


def template_table_data(data):
	result = ''
	for r in data:
		result += '<tr>'
		for d in r:
			result += '<td'
			result += ' class=\"' + d[1] + '\">' if d[1] != 'none' else '>'
			result += '<a href=\"'+ d[2] + '.html\">'
			result += str(d[0])
			result += '</a></td>'
			pass
		result += '</tr>'
	pass
	return result


def convert(names):
	result = []
	for name in names:
		result.append((name, 'none', ''))
	return result


def append_names(names, data):
	for i in range(len(data)):
		data[i].insert(0, names[i])
	pass
	return data


def append_x(names):
	names.insert(0, ('', 'none', ''))
	return names

def calcProtein(s1, s2, step):
	nseq1 = nseq2 = numbers = cursors = ''
	tmp = 0

	cnumber = step
	mismatch_count = 0

	for i in range(len(s1)):
		# print i
		if((i + 1) % step == 0):
			# print 'i'
			numbers += str(cnumber)
			tmp = len(str(cnumber)) - 1
			# print 'tmp', tmp

			cursors += '|'
			for x in range(len(str(cnumber)) - 1):
				# print '|', str(cnumber)
				cursors += ' '
			cnumber += step
		else:
			if(tmp != 0):
				tmp -= 1
			else:
				numbers += ' '
				cursors += ' '

		if(s1[i] != s2[i]):
			nseq1 += pystache.render('<span class=\"red\">{{symbol}}</span>', {'symbol': s1[i]})
			nseq2 += pystache.render('<span class=\"red\">{{symbol}}</span>', {'symbol': s2[i]})
			mismatch_count += 1
		else:
			nseq1 += s1[i]
			nseq2 += s2[i]

	return Bunch(nseq1=nseq1, nseq2=nseq2, numbers=numbers, cursors=cursors, mismatch_count=mismatch_count)

def calcNucleotide(s1, s2, step, mutations):
	nseq1 = nseq2 = numbers = cursors = ''
	tmp = 0

	cnumber = step
	mismatch_count = 0

	for i in range(len(s1)):
		# print i
		if((i + 1) % step == 0):
			# print 'i'
			numbers += str(cnumber)
			tmp = len(str(cnumber)) - 1
			# print 'tmp', tmp

			cursors += '|'
			for x in range(len(str(cnumber)) - 1):
				# print '|', str(cnumber)
				cursors += ' '
			cnumber += step
		else:
			if(tmp != 0):
				tmp -= 1
			else:
				numbers += ' '
				cursors += ' '

		if(s1[i] != s2[i]):
			color = 'yellow' if mutations[i] == 'nucleotide' else 'red'
			nseq1 += pystache.render('<span class=\"' + color + '\">{{symbol}}</span>', {'symbol': s1[i]})
			nseq2 += pystache.render('<span class=\"' + color + '\">{{symbol}}</span>', {'symbol': s2[i]})
			mismatch_count += 1
		else:
			nseq1 += s1[i]
			nseq2 += s2[i]

	return Bunch(nseq1=nseq1, nseq2=nseq2, numbers=numbers, cursors=cursors, mismatch_count=mismatch_count)

def createModelFotTemplate(seq1, seq2, alignment):
	dnaseq1 = alignment.aseq1
	dnaseq2 = alignment.aseq2

	aseq1 = alignment.aaseq1
	aseq2 = alignment.aaseq2

	dnalen1 = str(len(seq1.seq))
	dnalen2 = str(len(seq2.seq))

	sp1 = seq1.specie
	sp2 = seq2.specie

	fasta_headers_1 = [{'header':header} for header in seq1.row_headers]
	fasta_headers_2 = [{'header':header} for header in seq2.row_headers]

	b1 = calcNucleotide(dnaseq1, dnaseq2, 50, alignment.mutations)
	b2 = calcProtein(aseq1, aseq2, 20)

	return {
		'sp1': sp1,
		'sp2': sp2,
		'dnalen1': dnalen1,
		'dnalen2': dnalen2,
		'fasta_headers_1': fasta_headers_1,
		'fasta_headers_2': fasta_headers_2,
		'aseq1': b1.nseq1 + '     ',
		'aseq2': b1.nseq2 + '     ',
		'score': str(alignment.score),
		'mismatch_count': str(b1.mismatch_count),
		'numbers': b1.numbers + '     ',
		'cursors': b1.cursors + '     ',
		'aaseq1': b2.nseq1 + '     ',
		'aaseq2': b2.nseq2 + '     ',
		'amismatch_count': str(b2.mismatch_count),
		'anumbers': b2.numbers + '     ',
		'acursors': b2.cursors + '     ',
		'link': 'result.html'
	}

def copy_css():
	cssdir = './out/css' 
	if not os.path.exists(cssdir):
		os.makedirs(cssdir)

	shutil.copy2('./assets/style.css', './out/css')

def present(alignment_data):
	names = []
	sequences_count = len(alignment_data.sequences)

	copy_css()

	tmp = open('./assets/aligment_template.mustache', 'r')
	mustache = tmp.read()

	data = [[(0, 'none', '') for i in range(sequences_count)] for j in range(sequences_count)]
	for i in range(sequences_count):
		names.append(alignment_data.sequences[i].specie)

		for j in range(sequences_count):
			if (alignment_data.sequences[i].specie == alignment_data.sequences[j].specie):
				marker = 'green'
			elif (alignment_data.alignments[i][j].score < 8):
				marker = 'red'
			else:
				marker = 'none'

			out = open('out/' + str(i) + '-' + str(j) + '.html', 'w')
			html = pystache.render(mustache, createModelFotTemplate(alignment_data.sequences[i], alignment_data.sequences[j], alignment_data.alignments[i][j]))
			out.write(html)
			data[i][j] = (alignment_data.alignments[i][j].score, marker, str(i)+'-'+str(j))


	generate_template(data, names)

