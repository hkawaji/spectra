#!/bin/bash

fivePrimeFilterMinSigCount=0
fivePrimeFilterMinSigRatio=0.2
threePrimeFilterMinInternalPrimingRatio=0.5
matchRef=best
lessPrioritizeSingleCapCount=yes
intergenicRescue=yes
filterPath=./spectra/speFilter.sh
debug_dir=
nonSelectedFile=

usage ()
{
  cat <<EOF

  gunzip -c ISOGROUP_OUTPUT.bed.gz \\
  | $0 \\
       [ -n fileOfNotSelectedModels ]
       [ -c fivePrimeFilterMinSigCount ($fivePrimeFilterMinSigCount) ]
       [ -r fivePrimeFilterMinSigRatio ($fivePrimeFilterMinSigRatio) ]
       [ -i threePrimeFilterMinInternalPrimingRatio ($threePrimeFilterMinInternalPrimingRatio) ]
       [ -m matchRef (best[all|agnostic]) ]
       [ -l lessPrioritizeSingleCapCount (yes[no]) ]
       [ -u intergenicRescue (yes[no]) ]
       [ -f filterPath (${filterPath}) ]
       [ -d debug_dir ]

EOF
  exit 1
}

### handle options
while getopts n:c:r:i:m:l:u:f:d: opt
do
  case ${opt} in
  n) nonSelectedFile=${OPTARG};;
  c) fivePrimeFilterMinSigCount=${OPTARG};;
  r) fivePrimeFilterMinSigRatio=${OPTARG};;
  i) threePrimeFilterMinInternalPrimingRatio=${OPTARG};;
  m) matchRef=${OPTARG};;
  l) lessPrioritizeSingleCapCount=${OPTARG};;
  u) intergenicRescue=${OPTARG};;
  f) filterPath=${OPTARG};;
  d) debug_dir=${OPTARG};;
  *) usage;;
  esac
done

export LC_ALL=C

### setup tmpdir
tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15

awk '{print}' > ${tmpdir}/infile.bed

###
### Step 1
###
cat ${tmpdir}/infile.bed \
| $filterPath \
  -c $fivePrimeFilterMinSigCount \
  -r $fivePrimeFilterMinSigRatio \
  -i $threePrimeFilterMinInternalPrimingRatio \
  -m $matchRef \
| sort -k1,1 -k2,2n \
> ${tmpdir}/selected.bed

cat ${tmpdir}/infile.bed \
| $filterPath \
  -c $fivePrimeFilterMinSigCount \
  -r $fivePrimeFilterMinSigRatio \
  -i $threePrimeFilterMinInternalPrimingRatio \
  -m $matchRef \
  -v \
| sort -k1,1 -k2,2n \
> ${tmpdir}/discarded.bed



###
### Step 2 (less prioritize single evidence model when stronger ones are present)
###

if [ "$lessPrioritizeSingleCapCount" == "yes" ]; then

  cat ${tmpdir}/selected.bed \
  | awk 'BEGIN{OFS="\t"}{
      match($4, "capSigCount=[^,]+")
      if ( RLENGTH > 0) { capSigCount = substr($4, RSTART + 12, RLENGTH - 12) + 0 }
      print $0, capSigCount
    }' \
  > ${tmpdir}/selected.bed.tmp

  bedtools map -nonamecheck -f 0.5 -r -s -a ${tmpdir}/selected.bed.tmp -b ${tmpdir}/selected.bed.tmp -c 15 -o max \
  | awk --assign tmpdir=$tmpdir 'BEGIN{
      OFS="\t"
      discardFile = tmpdir "/discarded2.bed"
    }{
      match($4, "matchRef=[^,]+")
      if ( RLENGTH > 0) {
       matchRef = substr($4, RSTART + 9, RLENGTH - 9)
      }
      if ( ( $15 == 1 ) && ( $16 > 1 ) && (matchRef == "NA")) { print > discardFile }
      #if ( ( $15 == 1 ) && ( $16 > 1 ) ) { print > discardFile }
      else { print }
    }' \
  | cut -f 1-14 \
  > ${tmpdir}/selected2.bed

  cat \
    ${tmpdir}/discarded.bed \
    ${tmpdir}/discarded2.bed \
  | cut -f 1-14 \
  | sort -k1,1 -k2,2n \
  > ${tmpdir}/discarded2.bed.tmp
  mv -f ${tmpdir}/discarded2.bed.tmp ${tmpdir}/discarded2.bed

else

  cp -f ${tmpdir}/selected.bed ${tmpdir}/selected2.bed
  cp -f ${tmpdir}/discarded.bed ${tmpdir}/discarded2.bed

fi



###
### Step 3 (select gene models within uncovered region)
###

if [ "$intergenicRescue" == "yes" ]; then

  intersectBed -nonamecheck -wa -u -s \
    -a ${tmpdir}/discarded2.bed \
    -b ${tmpdir}/selected2.bed \
  > ${tmpdir}/discarded3.bed

  intersectBed -nonamecheck -wa -s -v \
    -a ${tmpdir}/discarded2.bed \
    -b ${tmpdir}/selected2.bed \
  > ${tmpdir}/discarded2_non_ovlp.bed

  bedtools map -nonamecheck -f 0.5 -r -s -a ${tmpdir}/discarded2_non_ovlp.bed -b ${tmpdir}/discarded2_non_ovlp.bed -c 5 -o max \
  | awk --assign tmpdir=$tmpdir 'BEGIN{
      OFS="\t"
      discardFile = tmpdir "/discarded3.bed"
    }{
      if( $5 == $15 ){ print } else { print >> discardFile }
    }' \
  | cut -f 1-14 \
  > ${tmpdir}/selected3.bed

  cat \
    ${tmpdir}/selected2.bed \
    ${tmpdir}/selected3.bed \
  | sort -k1,1 -k2,2n -k4,4 \
  > ${tmpdir}/selected3.bed.tmp
  mv -f ${tmpdir}/selected3.bed.tmp ${tmpdir}/selected3.bed

  cat ${tmpdir}/discarded3.bed \
  | cut -f 1-14 \
  | sort -k1,1 -k2,2n -k4,4  \
  > ${tmpdir}/discarded3.bed.tmp
  mv -f ${tmpdir}/discarded3.bed.tmp ${tmpdir}/discarded3.bed

else

  cp -f ${tmpdir}/selected2.bed ${tmpdir}/selected3.bed
  cp -f ${tmpdir}/discarded2.bed ${tmpdir}/discarded3.bed

fi




###
### output
###

cat ${tmpdir}/selected3.bed 

if [ "${nonSelectedFile}" != "" ]; then
  cat ${tmpdir}/discarded3.bed \
  > ${nonSelectedFile}
fi

if [ "${debug_dir-}" != "" ]; then
  mv ${tmpdir} ${debug_dir}
fi


