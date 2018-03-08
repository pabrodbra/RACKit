of = open("../contig.simple", 'r')
nf = open("../clean_contig.simple", 'w')

smallest_read = 200
current_content = []

for line in of:
	if len(current_content) == 2:
		if len(current_content[1]) >= smallest_read:
			nf.write(current_content[0])
			nf.write(current_content[1])
			current_content = []
	if len(current_content) < 2:
		current_content.append(line)
of.close()
nf.close()
