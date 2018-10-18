 ____   ____         _       ______        _       _______   
|_  _| |_  _|       / \    .' ____ \      / \     |_   __ \  
  \ \   / /______  / _ \   | (___ \_|    / _ \      | |__) | 
   \ \ / /|______|/ ___ \   _.____`.    / ___ \     |  ___/  
    \ ' /       _/ /   \ \_| \____) | _/ /   \ \_  _| |_     
     \_/       |____| |____|\______.'|____| |____||_____|    

DESCRIPTION
V-ASAP is a pipeline designed to assemble viral genome sequence from amplicon-seq reads.                                                              

DEPENDENCIES
	- Bowtie 2 version 2.2.6
	- cap3 VersionDate: 02/10/15
	- CASPER v0.8.2
	- CD-HIT version 4.6
	- cutadapt 1.16
	- FastQC v0.11.5
	- FASTX Toolkit 0.0.14
	- Python 3 & Biopython

INSTALLATION with Docker 
	1. install docker (https://docs.docker.com/install/linux/docker-ce/ubuntu/#docker-ee-customers)
	2. build the image : 
		$ cd /your_path/V-ASAP/
		$ sudo docker build --no-cache -t v-asap .
	3. run the container : 
		$ sudo docker run -it -v /path_of_your_dataset:/data v-asap /bin/bash 
	4. launch v-asap : 
		$ V-ASAP.sh -1 /data/reads1.fastq  -2 /data/reads2.fastq -p /data/primers.csv -o /data/V-ASAP_results

EXECUTION
$ V-ASAP.sh -1 R1.fastq  -2 R2.fastq -p primers.csv -o output_dir


