#!/bin/sh


export LC_ALL=C
SORT_OPT_BASE=" --compress-program=lzop --buffer-size=5G "

### setup tmpdir
tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15

bed12ToIntronCoords_geneSize_bed12 ()
{
  awk 'BEGIN{OFS="\t"}{
    chrom = $1; chromStart = $2; chromEnd = $3; 
    name = $4; score = $5; strand = $6; blockCount = $10;
    split($11,blockSizes,","); split($12,blockStarts,",");
    intronCoords = ""
    geneSize = chromEnd - chromStart
    for (i = 2; i<= blockCount; i++)
    {
      intronStart = chromStart + blockStarts[i-1] + blockSizes[i-1]
      intronEnd = chromStart + blockStarts[i]
      intronCoord = sprintf("%s:%s:%s:%s", chrom, intronStart, intronEnd, strand)
      intronCoords = intronCoords intronCoord ","
    }
    if (intronCoords != "") {
      print intronCoords, geneSize, $0
    }
  }'
}

selectTop ()
{
  awk 'BEGIN{OFS="\t";prevKey=""}{
    currKey = $1
    if ( prevKey != currKey ) {print}
    prevKey = currKey
  }'
}

usage ()
{
  cat <<EOF

usage: gunzip -c model.bedDetail.gz | $0 

EOF
  exit 1
}



#### handle options
while getopts i:h: opt
do
  case ${opt} in
  *) usage;;
  esac
done

awk --assign tmpdir=$tmpdir '{
    if($10 == 1){print > tmpdir "/res.bed12" }
    else{ print > tmpdir "/beforeSelect.bed12"}
  }'

cat ${tmpdir}/beforeSelect.bed12 \
| bed12ToIntronCoords_geneSize_bed12 \
| sort -k1,1 -k7,7nr -k2,2nr ${SORT_OPT_BASE} \
| cut -f 1,3- \
| selectTop \
| cut -f 2-  \
>> ${tmpdir}/res.bed12

sort -k1,1 -k2,2n ${SORT_OPT_BASE} ${tmpdir}/res.bed12


