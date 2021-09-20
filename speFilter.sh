#!/bin/bash

fivePrimeFilterMinSigCount=0
fivePrimeFilterMinSigRatio=0.2
threePrimeFilterMinInternalPrimingRatio=0.5
matchRef="best"
invert_match="false"

usage ()
{
  cat <<EOF

  gunzip -c ISOGROUP_OUTPUT.bed.gz \\
  | $0 \\
       [ -c fivePrimeFilterMinSigCount ($fivePrimeFilterMinSigCount) ]
       [ -r fivePrimeFilterMinSigRatio ($fivePrimeFilterMinSigRatio) ]
       [ -i threePrimeFilterMinInternalPrimingRatio ($threePrimeFilterMinInternalPrimingRatio) ]
       [ -m matchRef (best|all|agnostic) ]
       [ -v ] (for invert match, such as "grep -v")

EOF
  exit 1
}

### handle options
while getopts c:r:i:m:v opt
do
  case ${opt} in
  c) fivePrimeFilterMinSigCount=${OPTARG};;
  r) fivePrimeFilterMinSigRatio=${OPTARG};;
  i) threePrimeFilterMinInternalPrimingRatio=${OPTARG};;
  m) matchRef=${OPTARG};;
  v) invert_match="true";;
  *) usage;;
  esac
done

export LC_ALL=C

### setup tmpdir
tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15
touch ${tmpdir}/match.txt
touch ${tmpdir}/discard.bed



awk \
  --assign fivePrimeFilterMinSigCount=$fivePrimeFilterMinSigCount \
  --assign fivePrimeFilterMinSigRatio=$fivePrimeFilterMinSigRatio \
  --assign threePrimeFilterMinInternalPrimingRatio=$threePrimeFilterMinInternalPrimingRatio \
  --assign matchRef=$matchRef \
  --assign invert_match=$invert_match \
  --assign tmpdir=$tmpdir \
  'BEGIN{OFS="\t"; matchFile = tmpdir "/match.txt" }{

    match($4, "matchRef=[^,]+")
    if (( RLENGTH > 0) && ( matchRef != "agnostic")) {
      matchRef = substr($4, RSTART + 9, RLENGTH - 9)
      if ( matchRef != "NA" )
      {
        print matchRef , $0 > matchFile
        next
      }
    }

    match($4, "capSigCount=[^,]+")
    if ( RLENGTH > 0) { capSigCount = substr($4, RSTART + 12, RLENGTH - 12) + 0 }

    match($4, "capSigRatio=[^,]+")
    if ( RLENGTH > 0) { capSigRatio = substr($4, RSTART + 12, RLENGTH - 12) + 0 }

    match($4, "lastExonOverlapWithOtherInternalExons=[^,]+")
    if ( RLENGTH > 0) { lastExonOverlapWithOtherInternalExons = substr($4, RSTART + 38, RLENGTH - 38) + 0 }

    match($4, "internalPrimingRatio=[^,]+")
    if ( RLENGTH > 0) { internalPrimingRatio = substr($4, RSTART + 21, RLENGTH - 21) + 0 }

    flag = "yes"
    if ( capSigCount < fivePrimeFilterMinSigCount ) { flag = "no"}
    if ( capSigRatio < fivePrimeFilterMinSigRatio ) { flag = "no"}
    if ( ( internalPrimingRatio > threePrimeFilterMinInternalPrimingRatio ) &&
         ( lastExonOverlapWithOtherInternalExons > 0 ) ) { flag = "no"}

    # invert
    if ( invert_match == "true" ) {
      if ( flag == "yes" ) { flag = "no" } else { flag = "yes" }
    }

    # print
    if ( flag == "yes" ){ print } 
  }' \
> ${tmpdir}/out.bed


if [ "${matchRef}" == "all" ]; then

  if [ "${invert_match}" == "false" ]; then
    cat ${tmpdir}/match.txt \
    | cut -f 2- \
    >> ${tmpdir}/out.bed
  fi

elif [ "${matchRef}" == "best" ]; then

  cat ${tmpdir}/match.txt \
  | sort -k1,1 -k6,6nr \
  | awk --assign invert_match=$invert_match 'BEGIN{
      OFS="\t"
      prev_ref = ""
      prev_count = -1
    }{
      curr_ref = $1
      curr_count = $6
      flag = "no"
      if ( prev_ref != curr_ref ) { flag = "yes" }
      else if ( max_count == curr_count ) { flag = "yes"}

      # invert
      if ( invert_match == "true" ) {
        if ( flag == "yes" ) { flag = "no" } else { flag = "yes" }
      }

      # print
      if ( flag == "yes" ){ print }

      # set for next
      if ( prev_ref != curr_ref ) { max_count = curr_count }
      prev_ref = curr_ref
    }' \
  | cut -f 2- \
  >> ${tmpdir}/out.bed

fi


cat ${tmpdir}/out.bed \
| sort -k1,1 -k2,2n -k4,4 \
| uniq


