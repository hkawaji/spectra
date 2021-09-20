#!/bin/bash

fivePrimeFilterMinSigCount=0
fivePrimeFilterMinSigRatio=0.2
fivePrimeFilterMinSigCountIncl=0
fivePrimeFilterMinSigRatioIncl=0
threePrimeFilterMinInternalPrimingRatio=0.5

usage ()
{
  cat <<EOF

  gunzip -c ISOGROUP_OUTPUT.bed.gz \\
  | $0 \\
       [ -c fivePrimeFilterMinSigCount ($fivePrimeFilterMinSigCount) ]
       [ -r fivePrimeFilterMinSigRatio ($fivePrimeFilterMinSigRatio) ]
       [ -d fivePrimeFilterMinSigCountIncl ($fivePrimeFilterMinSigCountIncl) ]
       [ -s fivePrimeFilterMinSigRatioIncl ($fivePrimeFilterMinSigRatioIncl) ]
       [ -i threePrimeFilterMinInternalPrimingRatio ($threePrimeFilterMinInternalPrimingRatio) ]

EOF
  exit 1
}

### handle options
while getopts c:r:d:s:i: opt
do
  case ${opt} in
  c) fivePrimeFilterMinSigCount=${OPTARG};;
  r) fivePrimeFilterMinSigRatio=${OPTARG};;
  d) fivePrimeFilterMinSigCountIncl=${OPTARG};;
  s) fivePrimeFilterMinSigRatioIncl=${OPTARG};;
  i) threePrimeFilterMinInternalPrimingRatio=${OPTARG};;
  *) usage;;
  esac
done

awk \
  --assign fivePrimeFilterMinSigCount=$fivePrimeFilterMinSigCount \
  --assign fivePrimeFilterMinSigRatio=$fivePrimeFilterMinSigRatio \
  --assign fivePrimeFilterMinSigCountIncl=$fivePrimeFilterMinSigCountIncl \
  --assign fivePrimeFilterMinSigRatioIncl=$fivePrimeFilterMinSigRatioIncl \
  --assign threePrimeFilterMinInternalPrimingRatio=$threePrimeFilterMinInternalPrimingRatio \
  'BEGIN{OFS="\t"}{
    match($4, "matchRef=[^,]+")
    if ( RLENGTH > 0) { matchRef = substr($4, RSTART + 9, RLENGTH - 9)  }

    match($4, "capSigCount=[^,]+")
    if ( RLENGTH > 0) { capSigCount = substr($4, RSTART + 12, RLENGTH - 12) + 0 }

    match($4, "capSigRatio=[^,]+")
    if ( RLENGTH > 0) { capSigRatio = substr($4, RSTART + 12, RLENGTH - 12) + 0 }

    match($4, "capSigCountIncl=[^,]+")
    if ( RLENGTH > 0) { capSigCountIncl = substr($4, RSTART + 16, RLENGTH - 16) + 0 }

    match($4, "capSigRatioIncl=[^,]+")
    if ( RLENGTH > 0) { capSigRatioIncl = substr($4, RSTART + 16, RLENGTH - 16) + 0 }

    match($4, "lastExonOverlapWithOtherInternalExons=[^,]+")
    if ( RLENGTH > 0) { lastExonOverlapWithOtherInternalExons = substr($4, RSTART + 38, RLENGTH - 38) + 0 }

    match($4, "internalPrimingRatio=[^,]+")
    if ( RLENGTH > 0) { internalPrimingRatio = substr($4, RSTART + 21, RLENGTH - 21) + 0 }

    #print matchRef , capSigCount, capSigRatio, lastExonOverlapWithOtherInternalExons, internalPrimingRatio

    flag = "yes"
    if ( capSigCount < fivePrimeFilterMinSigCount ) { flag = "no"}
    if ( capSigRatio < fivePrimeFilterMinSigRatio ) { flag = "no"}
    if ( capSigCountIncl < fivePrimeFilterMinSigCountIncl ) { flag = "no"}
    if ( capSigRatioIncl < fivePrimeFilterMinSigRatioIncl ) { flag = "no"}
    if ( ( internalPrimingRatio > threePrimeFilterMinInternalPrimingRatio ) &&
         ( lastExonOverlapWithOtherInternalExons > 0 ) ) { flag = "no"}
    if ( flag == "yes" ){ print }
  }'
