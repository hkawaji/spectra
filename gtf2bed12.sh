#!/bin/bash

infile=$1

usage ()
{
  cat <<EOF
usage: $0 INFILE.gtf.gz

EOF
  exit 1
}

if [ ! -n "${infile-}" ]; then
  usage
fi


tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15

gunzip -c ${infile} > ${tmpdir}/tmp.gtf
gtfToGenePred ${tmpdir}/tmp.gtf ${tmpdir}/tmp.gp
genePredToBed ${tmpdir}/tmp.gp ${tmpdir}/tmp.bed
cat ${tmpdir}/tmp.bed

