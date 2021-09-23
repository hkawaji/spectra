#!/bin/bash

#infile=a.fa.gz
infile=$1

### setup for later
tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15

gunzip -c ${infile} \
| seqtk seq -F + - \
| gzip -c --fast \
> ${tmpdir}/infile.fq.gz

picard FastqToSam \
  -Xmx1g \
  --FASTQ  ${tmpdir}/infile.fq.gz  \
  --OUTPUT ${tmpdir}/infile.sam \
  --SORT_ORDER unsorted \
  --SAMPLE_NAME dummy
#rm -f ${tmpdir}/infile.fq.gz

python utils/flip_reads.py \
  --f ${tmpdir}/infile.sam \
  --o ${tmpdir}/infile.flip.sam
rm -f ${tmpdir}/infile.sam \

#cp -rp ${tmpdir} .

samtools fastq ${tmpdir}/infile.flip.sam \
| cutadapt --quiet -a "A{100}"  -o /dev/stdout - 


