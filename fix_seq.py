#!/usr/bin/env python
from Bio import SeqIO
import sys
for record in SeqIO.parse(sys.argv[1],"fasta"):
	print(">genecall_"+record.id+" ")
	print(record.seq+"*")
