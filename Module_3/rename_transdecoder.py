## Renames sequences from the TransDecoder output.

import sys

Fasta_file = sys.argv[1]


outfile=open("orciraptor_transdecoder.pep_renamed.fasta", "w")


Fasta_infile = open(Fasta_file)

lines = Fasta_infile.readlines()
for l in range(len(lines)-1):
	line = lines[l]
	#print compline
	if line[0] == ">":
		new_seq_name = ("%s_%s" % (line.split()[0].split("_")[6].strip(), line.split()[0].split("_")[7].strip()))
		outfile.write(">%s\n"%(new_seq_name))
		m = 1
		while not (lines[l+m][0] == ">"):
			if (l+m) == (len(lines)-1):
				line = lines[l+m]
				outfile.write(line)
				break
			else:
				next = lines[l+m]
				outfile.write(next)
				m = m+1


Fasta_infile.close()
outfile.close()
