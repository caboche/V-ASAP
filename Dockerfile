FROM ubuntu:18.04
MAINTAINER Florence MAURIER

RUN apt-get update ; apt-get upgrade -y ; apt-get install apt-utils -y ; apt-get install python3 python3-setuptools python3-biopython cutadapt bowtie2 fastx-toolkit fastqc cd-hit -y 

ADD dependencies/casper_v0.8.2/casper /usr/bin/.
ADD dependencies/CAP3/cap3 /usr/bin/.
ADD bin/print_most_represented_sequence_cd-hit.py /usr/bin/.
ADD bin/V-ASAP.sh /usr/bin/.
ADD bin/phiX174.fasta /usr/bin/.

