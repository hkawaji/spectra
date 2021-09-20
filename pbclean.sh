#!/bin/bash

usage ()
{
  cat <<EOF

  $0 -i infile(s) \\
     -s pb_same \\
     -a pb_adaptors \\
     [-f /path/to/flip_reads.py] \\
     [-p parallel(10)]


This script produce a cleaned-up fastq files from PacBio CCS raw reads
produced by a template switch-based cDNA library protocol [^1],
based on the steps described in its analysis protocol v3.1.

Dependencies:
* [flip_reads.py](https://github.com/mortazavilab/ENCODE-references)
* [isoseq3](https://github.com/PacificBiosciences/IsoSeq) tested in v3.4.0
* [lima](https://github.com/PacificBiosciences/barcoding) tested in v2.2.0
* [samtools](https://github.com/samtools/samtools) tested in v1.13


Motivation
----------
The analysis protocol v3.0 [^3] trims 5'-end too much,
while v3.1 [^2] does trim only common adaptor part, not GGG.
(probably due to 'isoseq3 refine'?, confimed in v3.4.0 and v3.2.2)
This script runs the v3.1 as it is in 'Part I', and trim
GGG stretch subsequently.

Note: the library protocol [^1] describes TSO sequence as below,
where 'common' part is also used as 3' adaptor too.

  5'-AAGCAGTGGTATCAACGCAGAGTACrGrG+G-3'
     |--- common ------------|

[^1]: Protocol to build non-size selected cDNA libraries for Pacific Biosciences
    long-read sequencing Version 3.0 (October, 2020)
    https://www.encodeproject.org/experiments/ENCSR507JOF/

[^2]: ENCODE Long Read RNA-Seq Analysis Protocol for Human Samples (v3.1)
  https://www.encodeproject.org/experiments/ENCSR634AKC/

[^3]: ENCODE Long Read RNA-Seq Analysis Protocol for Human Samples (v3.0)
  https://www.encodeproject.org/experiments/ENCSR507JOF/

EOF
  exit 1
}

pb_same=in/PB_adapters_same.fasta
pb_adapters=in/PB_adapters.fasta
parallel=10
flip_py=./utils/flip_reads.py

### handle options
while getopts i:p:s:a:f: opt
do
  case ${opt} in
  i) infiles=${OPTARG};;
  p) parallel=${OPTARG};;
  s) pb_same=${OPTARG};;
  a) pb_adapters=${OPTARG};;
  a) flip_py=${OPTARG};;
  *) usage;;
  esac
done
if [ ! -n "${infiles-}" ]; then usage; fi
if [ ! -n "${flip_py-}" ]; then usage; fi


### setup for later
tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15



part_I ()
{
  infile=$1
  prefix=$2

  ccs \
    -j ${parallel} \
    --noPolish \
    --minLength=10 \
    --minPasses=3 \
    --min-rq=0.9 \
    --min-snr=2.5 \
    --reportFile ${prefix}_ccs_report.txt \
    ${infile} \
    ${prefix}_ccs.bam

  lima ${prefix}_ccs.bam \
    ${pb_same} \
    ${prefix}_ccs_fl.bam \
    --ccs \
    --same \
    --num-threads ${parallel} \
    --min-score 80 \
    --min-end-score 0 \
    --min-signal-increase 10 \
    --min-score-lead 10 

  samtools view -h ${prefix}_ccs_fl.bam > ${prefix}_ccs_fl.sam
  python ${flip_py} --f ${prefix}_ccs_fl.sam --o ${prefix}_ccs_fl_flipped.sam
  samtools view -bS ${prefix}_ccs_fl_flipped.sam > ${prefix}_ccs_fl_flipped.bam

  isoseq3 refine \
    ${prefix}_ccs_fl_flipped.bam \
    ${pb_adapters} \
    ${prefix}_ccs_fl_flipped_clean.bam \
    --min-polya-length 20 \
    --require-polya \
    --num-threads ${parallel}

  samtools fastq ${prefix}_ccs_fl_flipped_clean.bam \
  | gzip -c \
  > ${prefix}_ccs_fl_flipped_clean.fq.gz
}


part_II ()
{
  infile=$1

  gunzip -c ${infile} \
  | grep --before 1 ^GGG  | grep ^@ | cut -b 2- \
  >  ${infile}.GGGlist

  seqtk subseq \
    ${infile} \
    ${infile}.GGGlist \
  | seqtk trimfq -b 3 - 
}

for infile in $infiles
do
  prefix=${tmpdir}/$(basename $infile .bam)
  part_I  ${infile} ${prefix}
  part_II ${prefix}_ccs_fl_flipped_clean.fq.gz
done
