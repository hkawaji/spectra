#!/bin/sh

SORT_OPT_BASE="--batch-size=100"

#chr1    1216913 1231935 m54284U_201027_092252/199160/ccs        60      -       1216913 1231935 255,0,0 7       775,176,159,114,137,479,44      0,1544,1855,6330,6918,11554,14978
bed12ToIntron ()
{
  awk 'BEGIN{OFS="\t"}{
    chrom = $1
    chromStart = $2
    chromEnd = $3
    name = $4
    score = $5
    strand = $6
    blockCount = $10
    split($11,blockSizes,",")
    split($12,blockStarts,",")
    for (i = 2; i<= blockCount; i++)
    {
      intronStart = chromStart + blockStarts[i-1] + blockSizes[i-1]
      intronEnd = chromStart + blockStarts[i]
      intronName = sprintf("%s:%s:%s:%s:%s,I%04d", name, chrom, chromStart, chromEnd, strand,i-1)
      print chrom, intronStart, intronEnd, intronName, 0, strand
    }
  }' 
}


most_freq_boundary ()
{
  awk 'BEGIN{OFS="\t"}{$4=$1":"$2":"$3":"$6;print}' \
  | mergeBed -s -c 4 -o collapse \
  | awk 'BEGIIN{OFS="\t"}{
    split("",counts)
    n=split($4,regions,",")
    for(i=1;i<=n;i++){ counts[ regions[i] ]++ }
    maxVal = 0
    maxKey = 0
    for(k in counts ) {if(counts[k] > maxVal){ maxVal = counts[k]; maxKey = k} }
    print maxKey ":" maxVal ":" n 
  }' \
  | awk 'BEGIN{FS=":"}{print $1"\t"$2"\t"$3"\t"$0"\t"$5"\t"$4}'
}

most_freq_intron ()
{
  local frac=$1
  local repeat=$2

  ### initial set
  cat ${tmpdir}/intron.bed \
  | most_freq_boundary \
  > ${tmpdir}/intron_mf.bed

  ### retrieve uncovered intron, due to overlap with others
  for rep in $(seq 1 ${repeat} )
  do
    intersectBed -v -s \
      -f $frac -r \
      -a ${tmpdir}/intron.bed \
      -b ${tmpdir}/intron_mf.bed \
    > ${tmpdir}/intron_uncovered.bed

    n=$( cat ${tmpdir}/intron_uncovered.bed | wc -l )
    if [[ $n -eq 0  ]]; then
      break
    fi

    cat ${tmpdir}/intron_uncovered.bed \
    | most_freq_boundary \
    > ${tmpdir}/intron_mf.bed.tmp

    cat ${tmpdir}/intron_mf.bed.tmp \
    >> ${tmpdir}/intron_mf.bed

    cat ${tmpdir}/intron_mf.bed \
    | sort -k1,1 -k2,2n $SORT_OPT_BASE \
    > ${tmpdir}/intron_mf.bed.tmp

    mv -f ${tmpdir}/intron_mf.bed.tmp ${tmpdir}/intron_mf.bed
  done

  if [[ $rep = $repeat ]]; then
    intersectBed -v -s \
      -f $frac -r \
      -a ${tmpdir}/intron.bed \
      -b ${tmpdir}/intron_mf.bed \
    > ${tmpdir}/intron_uncovered.bed
    n=$( cat ${tmpdir}/intron_uncovered.bed | wc -l )
    printf "[Warning] ${n} introns remain to be covered\n"  >&2
  fi
}



support_min_count=3
support_min_frac=0.95
mapQ=20
annotation_reference=
parallel=50
export LC_ALL=C

usage ()
{
  cat <<EOF

  $0 -i infile -g genome -o outprefix [-f support_min_frac ($support_min_frac)] [-c support_min_count ($support_min_count)] [-p parallel ($parallel)] [-a ANNOTATION_REFERENCE.bed12.gz]

EOF
  exit 1;
}


### handle options
while getopts i:g:o:c:f:q:a:p: opt
do
  case ${opt} in
  i) infile=${OPTARG};;
  g) genome=${OPTARG};;
  o) outprefix=${OPTARG};;
  c) support_min_count=${OPTARG};;
  f) support_min_frac=${OPTARG};;
  q) mapQ=${OPTARG};;
  a) annotation_reference=${OPTARG};;
  p) parallel=${OPTARG};;
  *) usage;;
  esac
done

if [ ! -n "${infile-}" ]; then usage; fi
if [ ! -n "${genome-}" ]; then usage; fi
if [ ! -n "${outprefix-}" ]; then usage; fi


### setup for later
tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15

minimap_opt_1st=
if [ ! -n "${annotation_reference-}" ]; then
  # do nothing
  :
else
  gunzip -c ${annotation_reference} \
  | bed12ToIntron \
  > ${tmpdir}/annRefIntron.bed
  minimap_opt_1st=" --junc-bed ${tmpdir}/annRefIntron.bed "
fi

###
### alignment, 1st round
###
minimap2 ${minimap_opt_1st} -t ${parallel} -a -x splice:hq ${genome} ${infile} \
| samtools view -b -q ${mapQ} -F 0x900 - \
> ${tmpdir}/first.bam

###
### intron prep
###
bamToBed -bed12 -i ${tmpdir}/first.bam \
| bed12ToIntron \
| sort -k1,1 -k2,2n $SORT_OPT_BASE \
> ${tmpdir}/intron.bed

most_freq_intron $support_min_frac 100

cat ${tmpdir}/intron_mf.bed \
| awk --assign support_min_count=$support_min_count \
  'BEGIN{OFS="\t"}{if ($5 >= support_min_count){print}}' \
> ${tmpdir}/intron_ref.bed


intersectBed -v -s \
  -f $support_min_frac -r \
  -a ${tmpdir}/intron.bed \
  -b ${tmpdir}/intron_ref.bed \
> ${tmpdir}/intron_unsupported.bed

###
### alignment, 2nd round
###
minimap2 \
  -t ${parallel} -a -x splice:hq \
  --junc-bed ${tmpdir}/intron_ref.bed \
  ${genome} ${infile} \
| samtools view -b -q ${mapQ} -F 0x900 - \
| samtools sort -  \
> ${outprefix}.bam
gzip -c ${tmpdir}/intron_ref.bed > ${outprefix}_intron_ref.bed.gz
#gzip -c ${tmpdir}/intron_unsupported.bed > ${outprefix}_intron_unsupported.bed.gz


