from string import Template

tmpl_file = open('assets/template.html', 'r')
output_file = open('result.html', 'w')

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
			result += str(d[0])
			result += '</td>'
			pass
		result += '</tr>'
	pass
	return result


def convert(names):
	result = []
	for name in names:
		result.append((name, 'none'))
	return result


def append_names(names, data):
	for i in xrange(len(data)):
		data[i].insert(0, names[i])
	pass
	return data


def append_x(names):
	names.insert(0, ('', 'none'))
	return names


def present(alignment_data):
	names = []
	sequences_count = len(alignment_data.sequences)

	data = [[(0, 'none') for i in xrange(sequences_count)] for j in xrange(sequences_count)]
	for i in xrange(sequences_count):
		names.append(alignment_data.sequences[i].specie)
		for j in xrange(sequences_count):
			if (alignment_data.sequences[i].specie == alignment_data.sequences[j].specie):
				marker = 'green'
			elif (alignment_data.alignments[i][j].score < 8):
				marker = 'red'
			else:
				marker = 'none'

			data[i][j] = (alignment_data.alignments[i][j].score, marker)

	generate_template(data, names)
