#!/bin/bash

###
### LICENSE
###


usage(){
echo "V-ASAP.sh -1 <reads1.fastq> -2 <reads2.fastq> -p <primers.tsv> [-o <output_directory> -t <number_of_threads> -m <memory(Mo)> -d <minimal_depth>]
	-1				fastq file of forward overlaping paired end reads [required]
	-2				fastq file of reverse overlaping paired end reads [required]
	-p	--primers	primers in tab-separated values format (with one header line and colomn 1 for amplicon number; column 2 for the forward primer sequence; column 3 for the reverse primer sequence) [required]
	-o	--output	output directory (will be created if not present) [default V-ASAP_results/]
	-t	--threads	number of threads to use for parallelisable steps [default 12]
	-m	--memory	limit of memory to use for CD-hit step [default none]
	-d	--depth		minimal number of reads by amplicon [default 100]"
}

while [ "$1" != "" ]; do
    case $1 in
        -1 )					shift
								INPUT_READS_1=$1
                                ;;
        -2 )					shift
								INPUT_READS_2=$1
                                ;;
        -p | --primers )		shift
								PRIMERS=$1
                                ;;
        -o | --output )			shift
								OUTPUT_DIR=$1
                                ;;
        -t | --threads )		shift
								NB_THREADS=$1
                                ;;
        -m | --memory )			shift
								MEMORY=$1
                                ;;
        -d | --depth )			shift
								MINIMAL_DEPTH=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

if [ "$NB_THREADS" == "" ]; then NB_THREADS=12 ; fi
if [ "$MEMORY" == "" ]; then MEMORY=65536 ; fi
if [ "$MINIMAL_DEPTH" == "" ]; then MINIMAL_DEPTH=100 ; fi
if [ "$OUTPUT_DIR" == "" ]; then OUTPUT_DIR="V-ASAP_results" ; fi
if [ ! -d "$OUTPUT_DIR" ]; then mkdir $OUTPUT_DIR ; fi

PRIMER_ID_LIST=`tail -n +2 $PRIMERS | cut -f 1`

### MERGING READS ###
MERGED=${OUTPUT_DIR}/merged_reads
casper $INPUT_READS_1 $INPUT_READS_2 -l -o $MERGED -j -t ${NB_THREADS} >& ${OUTPUT_DIR}/log_merging_casper.txt

### CLEANING READS (ie PhiX sequence filtering) ###
PHIX_FASTA=`dirname $0`/phiX174.fasta
PHIX_INDEX=`dirname $0`/phiX174
if [ ! -f "${PHIX_INDEX}.1.bt2" ] ; then bowtie2-build ${PHIX_FASTA} ${PHIX_INDEX} ; fi 
INPUT_READS=${OUTPUT_DIR}/clean_reads.fastq
bowtie2 --local -p ${NB_THREADS} -x ${PHIX_INDEX} -S tmp.sam -U ${MERGED}.fastq --un ${INPUT_READS} >& ${OUTPUT_DIR}/log_clean_phiX.txt
rm tmp.sam

### GROUP READS BY AMPLICON ###
READS_DIR="${OUTPUT_DIR}/reads_by_amplicon/"
if [[ ! -d $READS_DIR ]] ; then  mkdir $READS_DIR ; fi 
for i in $PRIMER_ID_LIST
		do FORWARD=`grep ^${i}$'\t' $PRIMERS | cut -f 2`
		REVERSE=`grep ^${i}$'\t' $PRIMERS | cut -f 3`
		F_THRESHOLD=$((`echo $FORWARD | wc -c` - 1 ))
		R_THRESHOLD=$((`echo $REVERSE | wc -c` - 1 ))
		FORWARD_RC=`echo $FORWARD | tr "[ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv]" "[TACGAtacgaNnRrYySsWwMmKkVvHhDdBb]" | rev`
		REVERSE_RC=`echo $REVERSE | tr "[ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv]" "[TACGAtacgaNnRrYySsWwMmKkVvHhDdBb]" | rev`
		AMP_READS_TRIM=${READS_DIR}/"reads_${i}.fastq"
		cat $INPUT_READS | cutadapt -j ${NB_THREADS} --discard-untrimmed -g $FORWARD -a $FORWARD_RC -O $F_THRESHOLD - | cutadapt -j ${NB_THREADS} --discard-untrimmed -g $REVERSE -a $REVERSE_RC -O $R_THRESHOLD - > $AMP_READS_TRIM
done

### CLUSTERING STEP and SELECTION OF REPRESENTATIVE READS ###
CLUSTERS_DIR="${OUTPUT_DIR}/clusters-cdhit/"
if [[ ! -d ${CLUSTERS_DIR} ]] ; then  mkdir ${CLUSTERS_DIR} ; fi 
ONE_READ_PER_AMPLICON="${CLUSTERS_DIR}/one-read-per-amplicon.fastq"
if [[ -f "${ONE_READ_PER_AMPLICON}" ]] ;  then rm ${ONE_READ_PER_AMPLICON} ; fi 
for i in $PRIMER_ID_LIST
	do AMP_READS=${READS_DIR}/"reads_${i}.fastq"
	if  [ `grep ^+$ ${AMP_READS} | wc -l` -gt $MINIMAL_DEPTH ]
		then CLUSTER_CDHIT="${CLUSTERS_DIR}/clusters_cd-hit_${i}.fastq"
		cdhit-est -c 1 -d 0 -s 1 -i ${AMP_READS} -o ${CLUSTER_CDHIT} -M $MEMORY -T ${NB_THREADS} >& ${CLUSTERS_DIR}/log_cd-hit_${i}.txt
		print_most_represented_sequence_cd-hit.py ${CLUSTER_CDHIT} 1 | sed s/^"@".*/"@"${i}/g >> ${ONE_READ_PER_AMPLICON}
	fi 		
done

### ASSEMBLY STEP ###
ASS_DIR="${OUTPUT_DIR}/assembly-CAP3/"
if [[ ! -d $ASS_DIR ]] ; then  mkdir $ASS_DIR ; fi 
FASTQ="${OUTPUT_DIR}/clusters-cdhit/one-read-per-amplicon.fastq"
FASTA="${ASS_DIR}"/reads.fasta
fastq_to_fasta -n -i $FASTQ -o $FASTA
cap3 $FASTA -o 16 -s 251 -j 31 -p 66 -i 21  > ${ASS_DIR}/assembly-CAP3.aln
cat ${FASTA}.cap.contigs ${FASTA}.cap.singlets > `echo $FASTA | sed s/reads/assembly/g`
