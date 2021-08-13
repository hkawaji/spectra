#!/bin/bash


export LC_ALL=C
SORT_OPT_BASE="--batch-size=100"


### setup tmpdir
tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15



bamAddDownstreamNucAsSeqNameSuffix ()
{
  local infile=$1
  local genome=$2
  local targetLen=$3

  bamToBed -bed12 -i ${infile} \
  | awk --assign targetLen=$targetLen 'BEGIN{OFS="\t"}{
      chrom = $1
      chromStart = $2
      chromEnd = $3
      name = $4
      strand = $6
      if ( strand == "+"){
        chromStart = chromEnd
        chromEnd = chromEnd + targetLen
      } else {
        chromEnd = chromStart
        chromStart = chromEnd - targetLen
        if (chromStart < 0) {chromStart = 0}
      }
      print chrom, chromStart, chromEnd, name, 0, strand
    }' \
  | bedtools getfasta -tab -nameOnly -s -fi ${genome} -bed stdin \
  | awk 'BEGIN{OFS="\t"}{ $1 = substr($1,0,length($1) - 3); print }' \
  | sort -k1,1 $SORT_OPT_BASE \
  > ${tmpdir}/downNuc.txt

  samtools view -H ${infile} > ${tmpdir}/downNuc.sam
  samtools view ${infile} \
  | sort -k1,1 $SORT_OPT_BASE \
  | join -t "	" - ${tmpdir}/downNuc.txt \
  | awk 'BEGIN{OFS="\t"}{$1 = $1 ".3endDown-" $NF; NF = NF - 1; print}' \
  >> ${tmpdir}/downNuc.sam

  samtools view -b ${tmpdir}/downNuc.sam \
  | samtools sort - 
}


bamAddPolyASignalAsSeqNameSuffix ()
{
  # Nucleic Acids REs, 33:201-212, 2005
  # Top 5 PAS hexamers AAUAAA|AUUAAA|UAUAAA|AGUAAA|AAGAAA
  # Top 2 PAS hexamers AAUAAA|AUUAAA
  samtools view -h - \
  | awk 'BEGIN{OFS="\t"; terminalLen=50}
    function revcomp(seq)
    {
      a["T"]="A";a["A"]="T";a["C"]="G";a["G"]="C";a["N"]="N";
      a["t"]="A";a["a"]="T";a["c"]="G";a["g"]="C";a["n"]="N";
      seq_rc = ""
      for(i=length(seq);i!=0;i--) {
        seq_rc = seq_rc a[ substr(seq,i,1) ]
      }
      return seq_rc
    }
    {
      softClipLen = 0
      if ( match($0,/^@/) ) # header
      {
        # do nothing

      } else if ( and( $2 , 16 ) ) { # reverse
        match($6,/^([0-9]+)S/,buf)
        if (RLENGTH >= 0) { softClipLen = buf[1] }
        terminalSeq = substr($10, softClipLen, terminalLen)
        terminalSeq = revcomp( terminalSeq )
        match(terminalSeq,/(AATAAA|ATTAAA)/,buf)
        if (RLENGTH >= 0) { $1 = $1 ".3endPas" (RSTART - terminalLen) "-" buf[1]}
        else { $1 = $1 ".3endPas--NOTFOUND"} 
      } else {                       # forward
        match($6,/([0-9]+)S$/,buf)
        if (RLENGTH >= 0) { softClipLen = buf[1] }
        terminalSeq = substr($10,length($10) - softClipLen - terminalLen, terminalLen)
        match(terminalSeq,/(AATAAA|ATTAAA)/,buf)
        if (RLENGTH >= 0) { $1 = $1 ".3endPas" (RSTART - terminalLen) "-" buf[1] }
        else { $1 = $1 ".3endPas--NOTFOUND"} 
      }
      if (softClipLen <= 10) {
        print $0
      } else {
        printf "ReadFilteredOut: softClipLen %d\t%s\n", softClipLen, $0 > "/dev/stderr"
      }
    }' \
  | samtools view -b -
}
   


bamAdd5endAsSeqNameSuffix ()
{
  samtools view -h - \
  | awk 'BEGIN{OFS="\t"}
    function revcomp(seq)
    {
      a["T"]="A";a["A"]="T";a["C"]="G";a["G"]="C";a["N"]="N";
      a["t"]="A";a["a"]="T";a["c"]="G";a["g"]="C";a["n"]="N";
      seq_rc = ""
      for(i=length(seq);i!=0;i--) {
        seq_rc = seq_rc a[ substr(seq,i,1) ]
      }
      return seq_rc
    }
    {
      softClipLen = 0
      if ( match($0,/^@/) ) # header
      {
        # do nothing

      } else if ( and( $2 , 16 ) ) { # reverse
        match($6,/([0-9]+)S$/,buf)
        if (RLENGTH>=0){
          unmatchN = substr($10, length($10) - buf[1] + 1) 
          matchN = substr($10,length($10) - buf[1], 1)
          str = "5end-" revcomp(unmatchN) "-" revcomp(matchN)
          softClipLen = buf[1]
        } else {
          matchN = substr($10,length($10) )
          str = "5end--" revcomp(matchN)
        }
        $1 = $1 "." str

      } else {                       # forward
        match($6,/^([0-9]+)S/,buf)
        if (RLENGTH>=0){
          unmatchN = substr($10,0,buf[1])
          matchN = substr($10,buf[1]+1,1)
          str = "5end-" unmatchN "-" matchN
          softClipLen = buf[1]
        } else {
          matchN = substr($10,0,1)
          str = "5end--" matchN
        }
        $1 = $1 "." str
      }
      if (softClipLen <= 10) {
        print $0
      } else {
        printf "ReadFilteredOut: softClipLen %d\t%s\n", softClipLen, $0 > "/dev/stderr"
      }
    }' \
  | samtools view -b -
}


bamThreePrimeFilter ()
{
  local minRatioA=$1

  samtools view -h - \
  | awk --assign minRatioA=$minRatioA \
    'BEGIN{OFS="\t"}
    {
      if ( match($0,/^@/) ) # header
      {
        # do nothing
        print $0
      } else {
        match($1,/\.3endDown-([A-Za-z]*)/,buf)
        downNuc = buf[1]
        match($1,/\.3endPas-([0-9]*)-([A-Za-z]*)/,buf)
        pas = buf[2]
        countA = gsub("A|a","=",downNuc)
        ratioA = countA / ( length(downNuc) )
        if ( ( pas == "NOTFOUND") && (ratioA >= minRatioA) ) {
          printf "ReadFilteredOut: internalPriming\t%s\n", $0 > "/dev/stderr"
        } else { print $0 }
      }
    }' \
  | samtools view -b -
}


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
      intronName = sprintf("%s:%s:%s:%s:%s,I%04d:%04d", name, chrom, chromStart, chromEnd, strand, i-1, blockCount - 1)
      print chrom, intronStart, intronEnd, intronName, 0, strand
    }
  }' 
}



most_freq_intron_sub()
{
  local frac=$1
  local bedFreq=$( mktemp --tmpdir=${tmpdir} )

  awk 'BEGIN{OFS="\t"}{
    print $1":"$2":"$3":"$6
  }' \
  | sort ${SORT_OPT_BASE} \
  | uniq -c  \
  | awk '{print $2":"$1}' \
  | awk 'BEGIN{FS=":";OFS="\t"}{
      name=$1":"$2":"$3":"$4
      print $1,$2,$3,name,$5,$4
    }' \
  | sort -k1,1 -k2,2n ${SORT_OPT_BASE} \
  > ${bedFreq}

  bedtools map -s -f $frac -r \
    -a ${bedFreq} -b ${bedFreq} \
    -c 5,5,4 -o max,collapse,collapse \
  | awk '{if($5 == $7){print}}' \
  | awk 'BEGIN{OFS="\t"}{
      chrom = $1; chromStart = $2; chromEnd = $3;
      name = $4; score = $5; strand = $6;
      scoresN = split($8,scores,",")
      namesN  = split($9,names,",")
      flag = "yes"
      for (i=1;i<=scoresN;i++)
      {
        if ( scores[i] == score ) {
          split( names[i] , buf, ":")
          if ( buf[2] < chromStart ) {
            flag = "no"
          } else if (( buf[2] == chromStart ) && ( buf[3] < chromEnd )) {
            flag = "no"
          }
        }
      }
      if (flag == "yes") { print }
    }'
}


most_freq_intron()
{
  local infile=$1
  local frac=$2
  local bedMf=$( mktemp --tmpdir=${tmpdir} )

  cat ${infile} \
  | most_freq_intron_sub $frac \
  > ${bedMf}

  intersectBed -v -s -f $frac -r -a ${infile} -b ${bedMf} \
  > ${bedMf}.uncovered

  local n=$( cat ${bedMf}.uncovered | wc -l )
  if [[ $n -ne 0  ]]; then
    cat ${bedMf}.uncovered \
    | most_freq_intron_sub $frac \
    >> ${bedMf}
  fi

  cut -f 1-6 ${bedMf} \
  | awk 'BEGIN{OFS="\t"}{$4 = $4 ":" $5;print}'
}



intron_readsTobed12 ()
{
  awk 'BEGIN{OFS="\t"}{
    name = $2
    readsN = split(name,reads,",")

    #
    # define chromStart/chromEnd
    #
    split("",chromStartCounts,",")
    split("",chromEndCounts,",")
    for (i=1; i <= readsN; i++)
    {
      bufN = split(reads[i],buf,":")
      _chromStart = buf[3]
      _chromEnd = buf[4]
      chromStartCounts[ _chromStart ]++
      chromEndCounts[ _chromEnd ]++
    }

    chromStartMostFreq = 10000000000
    chromEndMostFreq = -1
    chromStartMostLeft = 10000000000
    chromEndMostRight = -1
    chromStartCountMax=0
    chromEndCountMax=0
    for (i=1; i <= readsN; i++)
    {
      bufN = split(reads[i],buf,":")
      readName = buf[1]
      chrom = buf[2]
      _chromStart = buf[3]
      _chromEnd = buf[4]
      strand = buf[5]

      if ( ( chromStartCounts[ _chromStart] >= chromStartCountMax ) &&
           ( _chromStart < chromStartMostFreq ) ) {
        chromStartMostFreq = _chromStart
        chromStartCountMax = chromStartCounts[ _chromStart]
      }
      if ( _chromStart < chromStartMostLeft ) { chromStartMostLeft = _chromStart }

      if ( ( chromEndCounts[ _chromEnd] >= chromEndCountMax ) &&
           ( _chromEnd > chromEndMostFreq ) ) {
        chromEndMostFreq = _chromEnd
        chromEndCountMax = chromEndCounts[ _chromEnd]
      }
      if ( _chromEnd > chromEndMostRight ) { chromEndMostRight = _chromEnd }
    }
    chromStart = chromStartMostFreq
    chromEnd = chromEndMostFreq

    # define blockSizes/blockStarts
    # chr10:100354632:100356531:+:2:2, ... ,chr10:100356764:100360733:+:2:2
    intronsN = split($1, introns, ",")
    blockCount = intronsN + 1
    blockSizes = ""
    blockStarts = "0,"

    # --- exceptional cases ---
    bufN = split(introns[1],buf,":")
    if ( buf[2] < chromStart ) { chromStart = buf[2] - 1}
    # -------------------------

    _end = chromStart
    for (i=1; i<= intronsN; i++)
    {
      bufN = split(introns[i],buf,":")
      _prevEnd   = _end
      _start = buf[2]
      _end = buf[3]
      blockSizes = blockSizes ( _start - _prevEnd ) ","
      blockStarts = blockStarts ( _end - chromStart ) ","
    }

    # --- exceptional cases ---
    if ( chromEnd < _end ) { chromEnd = _end + 1}
    # -------------------------

    blockSizes = blockSizes ( chromEnd - _end ) ","
    print chrom, chromStart, chromEnd, name, readsN, strand, 
      chromStart, chromStart, "0,0,0",
      blockCount, blockSizes, blockStarts
  }'
}


bed12ToBed12detail ()
{
  awk 'BEGIN{OFS="\t"}{
    buf = $4
    $4 = sprintf("SM%08d",NR)
    print $0, $4, buf
  }'
}



# [localNames]
# I0001:0012,I0002:0012, ... ,I0012:0012
# [refNames]
# chrX:154349036:154349361:-:1:1,chrX:154349565:154349648:-:1:1, ... ,chrX:154353457:154353553:-:1:1
intronSetFilter ()
{
  local minRefCount=$1
  awk --assign minRefCount=$minRefCount 'BEGIN{OFS="\t"}{
    localNamesN=split($2,localNames,",")
    refNamesN=split($3,refNames,",")
    flag = "ok"
    for (i=1;i<=localNamesN;i++){
      str = sprintf("I%04d:%04d",i,localNamesN)
      if ( localNames[i] != str ) { flag = "FileredOut:" localNames[i] " != " str}
    }
    for (i=1;i<=refNamesN;i++){
      split( refNames[i] , buf, ":")
      refCount = buf[5]
      if ( refCount < minRefCount ) {flag = "ReadFilteredOut: count of " refNames[i] " is < " minRefCount}
    }
    if (flag != "ok") {
      printf "%s\t%s\n", flag, $0 > "/dev/stderr"
    } else {
      print
    }
  }'
}

# input - read_intronLocalNames_intronRefNames.txt
# m54284U_200720_151958/100008005/ccs...:chr4:73404289:73421253:+   ...  ...
# output - read_intronLocalNames_intronRefNames_boundayMf.txt
# m54284U_200720_151958/100008005/ccs...:chr4:73404289:73421253:+   ...  ...  chr4:000:111:+

most_freq_boundary () {
  # read_intronLocalNames_intronRefNames.txt
  local infile=$1
  local frac=$2
  cat $infile \
  | awk 'BEGIN{OFS="\t"}{
      split($1,buf,":")
      chrom=buf[2];chromStart=buf[3];chromEnd=buf[4];strand=buf[5];
      print chrom, chromStart, chromEnd, $1, 0, strand
    }' \
  | sort -k1,1 -k2,2n \
  > ${infile}.tmpboundary

  bedtools map -s -f $frac -r \
    -a ${infile}.tmpboundary -b ${infile}.tmpboundary \
    -c 2,3 -o mode,mode \
  | awk 'BEGIN{OFS="\t"}{
      boundaryMf = $1 ":" $7 ":" $8 ":" $6
      print $4, boundaryMf
    }'\
  | sort -k1,1 \
  > ${infile}.tmpboundaryMf

  join $infile ${infile}.tmpboundaryMf \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}'
}




fivePrimeFilter ()
{
  local minSigCount=$1
  local minSigRatio=$2
  awk \
    --assign minSigCount=$minSigCount \
    --assign minSigRatio=$minSigRatio \
  'function isSig(seq) {
    if ( match(seq,/^([GCgc]*[Gg])$/,buf) )
    {
      gN = gsub("G|g", "=", buf[1])
      cN = gsub("C|c", "=", buf[1])
      if ( gN >= cN ) {return "yes"} else {return "no"}
    } else { return "no" }

  } BEGIN{ OFS="\t" } {
    readsN=split($14,reads,",")
    sigCount = 0
    for (i=1;i<=readsN;i++)
    {
      match(reads[i],/\.5end-([ATGCN]*)-([ATGCN]*)/,buf)
      unmatch = buf[1]
      alignmentStart = buf[2]
      if ( isSig(unmatch) == "yes" ) {sigCount++}
    }
    sigRatio = sigCount / readsN
    $4 = sprintf("%s,capSigCount=%d",$4,sigCount)
    $4 = sprintf("%s,capSigRatio=%.2f",$4,sigRatio)
    if ( ( sigCount >= minSigCount ) || ( sigRatio >= minSigRatio ) )
    {
      print
    }else{
      printf "ModelFilteredOut: capSigRatio %f, sigCount %f\t%s\n", sigRatio, sigCount, $0 > "/dev/stderr"
    }
  }'
}


#m54284U_200720_151958/101714429/ccs.3endDown-AAAAAAAAAAAAAAATTAGC.3endPas--NOTFOUND.5end-CG-G:chr5:138002847:138019921:-
threePrimeFilter ()
{
  local minRatioA=$1
  local minInternalPrimingRatio=$2

  awk \
    --assign minRatioA=$minRatioA \
    --assign minInternalPrimingRatio=$minInternalPrimingRatio \
  'BEGIN{OFS="\t"}{
    readsN=split($14,reads,",")
    internalPriming = 0
    for (i=1;i<=readsN;i++) {
      match(reads[i],/\.3endDown-([A-Za-z]*)/,buf)
      downNuc = buf[1]
      match(reads[i],/\.3endPas-([0-9]*)-([A-Za-z]*)/,buf)
      pas = buf[2]

      countA = gsub("A|a","=",downNuc)
      ratioA = countA / ( length(downNuc) )
      if ( ( pas == "NOTFOUND") && (ratioA >= minRatioA) ) { internalPriming ++ }
    }

    internalPrimingRatio = internalPriming / readsN
    $4 = sprintf("%s,internalPrimingRatio=%.2f",$4,internalPrimingRatio)
    if ( internalPrimingRatio >= minInternalPrimingRatio ) 
    {
      printf "ModelFilteredOut: internal priming ratio %f\t%s\n", internalPrimingRatio, $0 > "/dev/stderr"
    } else {
      print
    }
  }'
}



###
### main
###

mapQ=20
support_min_frac_intron=0.99
support_min_frac_boundary=0.9

usage ()
{
  cat <<EOF

  $0 -i infile -g genome [-q mapQ(${mapQ})] [-f support_min_frac_intron(${support_min_frac_intron})] [-r support_min_frac_boundary(${support_min_frac_boundary})]

EOF
  exit 1
}

### handle options
while getopts i:g:q:f:r: opt
do
  case ${opt} in
  i) infile=${OPTARG};;
  g) genome=${OPTARG};;
  q) mapQ=${OPTARG};;
  f) support_min_frac_intron=${OPTARG};;
  r) support_min_frac_boundary=${OPTARG};;
  *) usage;;
  esac
done

if [ ! -n "${infile-}" ]; then usage; fi
if [ ! -n "${genome-}" ]; then usage; fi

samtools view -bq ${mapQ} ${infile} > ${tmpdir}/infile.bam

bamAddDownstreamNucAsSeqNameSuffix ${tmpdir}/infile.bam ${genome} 20 \
| bamAddPolyASignalAsSeqNameSuffix \
| bamAdd5endAsSeqNameSuffix \
| bamToBed -bed12 -i stdin  \
| bed12ToIntron \
| sort -k1,1 -k2,2n $SORT_OPT_BASE \
> ${tmpdir}/intron.bed

most_freq_intron ${tmpdir}/intron.bed $support_min_frac_intron \
| sort -k1,1 -k2,2n $SORT_OPT_BASE \
> ${tmpdir}/intron_mf.bed

intersectBed -sorted -s -wa -wb \
  -f $support_min_frac_intron -r \
  -a ${tmpdir}/intron.bed \
  -b ${tmpdir}/intron_mf.bed \
> ${tmpdir}/intron_assignment.bed \

cat ${tmpdir}/intron_assignment.bed \
| sed -e 's/,/\t/' \
| sort -k4,4 -k5,5 $SORT_OPT_BASE \
| groupBy -g 4 -c 5,11 -o collapse \
> ${tmpdir}/read_intronLocalNames_intronRefNames.txt

most_freq_boundary \
  ${tmpdir}/read_intronLocalNames_intronRefNames.txt \
  $support_min_frac_boundary \
> ${tmpdir}/read_intronLocalNames_intronRefNames_boundaryMf.txt

cat ${tmpdir}/read_intronLocalNames_intronRefNames_boundaryMf.txt \
| intronSetFilter 1 \
| sort -k3,3 -k4,4 $SORT_OPT_BASE \
| groupBy -g 3,4 -c 1 -o collapse \
> ${tmpdir}/introns_boundaryMf_reads.txt

cat ${tmpdir}/introns_boundaryMf_reads.txt \
| cut -f 1,3 \
| intron_readsTobed12 \
| bed12ToBed12detail \
| fivePrimeFilter 3 0.5 \
| threePrimeFilter 0.5 0.5

