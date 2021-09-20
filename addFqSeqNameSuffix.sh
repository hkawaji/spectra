#!/bin/bash

suffix=$1

if [ ! -n "${suffix-}" ]; then
  cat <<EOF

usage: cat FASTQ | $0 suffix

EOF
  exit 1
fi

awk --assign suffix=$suffix '{
    if ((NR % 4) == 1) {
      $0 = $0 suffix
    }
    print
  }'


