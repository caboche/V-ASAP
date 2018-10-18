#!/usr/bin/python
# $ print_most_represented_sequence_cd-hit.py <clusters_cd-hit.fastq>
# (the file <clusters_cd-hit.fastq.clstr> should be in the same directory as <clusters_cd-hit.fastq>) 
import sys
import re
import os

def parse_clstr_file(clstr_adr):
	clstr = {}
	clstr_id = None
	clstr_file = open(clstr_adr, "r")
	line = clstr_file.readline()
	while line and line!="":
		if re.match("^>Cluster", line):
			clstr_id = line.split(" ")[1]
			clstr[clstr_id] = {}
			clstr[clstr_id]["nb"] = 0
			line = clstr_file.readline()
			while line and len(line)>0 and line[0]!=">":
				clstr[clstr_id]["nb"] += 1
				if re.match(".*\*$", line):
					clstr[clstr_id]["seq"] = line.split(">")[1][:-6]
				line = clstr_file.readline()
	return clstr

def find_most_representative_sequence(fastq_adr, clstr):
	max_clstr_size = 0
	max_clstr_id = None
	for clstr_id in clstr.keys():
		if clstr[clstr_id]["nb"] > max_clstr_size:
			max_clstr_size =  clstr[clstr_id]["nb"]
			max_clstr_id = clstr_id
	fastq_file = fastq_adr.split("/")
	sys.stderr.write(str(max_clstr_size)+"read(s) for file "+fastq_file[len(fastq_file)-1]+"\n")
	return clstr[max_clstr_id]["seq"]

def main():
	fastq_adr = sys.argv[1]
	clstr_adr = fastq_adr+".clstr"
	clstr = parse_clstr_file(clstr_adr)
	if (len(clstr) > 0) :
		seq = find_most_representative_sequence(fastq_adr, clstr)
		cmd = "sed s/\" .*$\"//g "+fastq_adr+" | grep \"^@"+seq+"$\" -A 3 "
		os.system(cmd)

if __name__ == "__main__":
    main()

