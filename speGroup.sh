#!/bin/bash


export LC_ALL=C
#SORT_OPT_BASE="--batch-size=100"
SORT_OPT_BASE=

### setup tmpdir
tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15


bamAddDownstreamNucAsSeqNameSuffix ()
{
  local infile=$1
  local genome=$2
  local targetLen=$3

  samtools view ${infile} \
  | cut -f 1 \
  | sort -k1,1 $SORT_OPT_BASE \
  > ${tmpdir}/seqids.txt

  bamToBed -bed12 -i ${infile} \
  | awk --assign targetLen=$targetLen 'BEGIN{OFS = "\t"}
    {
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
  | join -t "	"  -a 1  -e N -o 1.1,2.2 ${tmpdir}/seqids.txt - \
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
  | awk --assign tmpdir=${tmpdir} 'BEGIN{OFS="\t"; terminalLen=50}
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
        stde = tmpdir "/err.bamAddPolyASignalAsSeqNameSuffix.txt"
        printf "ReadFilteredOut: softClipLen %d\t%s\n", softClipLen, $0 >> stde
      }
    }' \
  | samtools view -b -
}
   


bamAdd5endAsSeqNameSuffix ()
{
  samtools view -h - \
  | awk --assign tmpdir=${tmpdir} 'BEGIN{OFS="\t"}
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
        stde = tmpdir "/err.bamAdd5endAsSeqNameSuffix.txt"
        printf "ReadFilteredOut: softClipLen %d\t%s\n", softClipLen, $0 >> stde
      }
    }' \
  | samtools view -b -
}


bamThreePrimeFilter ()
{
  local minRatioA=$1

  samtools view -h - \
  | awk --assign minRatioA=$minRatioA --assign tmpdir=${tmpdir} \
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
          stde = tmpdir "/err.bamThreePrimeFilter.txt"
          printf "ReadFilteredOut: internalPriming\t%s\n", $0 >> stde
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
  awk --assign tmpdir=${tmpdir} 'BEGIN{OFS="\t"}{
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

    outflag = "ok"
    _end = chromStart
    for (i=1; i<= intronsN; i++)
    {
      bufN = split(introns[i],buf,":")
      _prevEnd   = _end
      _start = buf[2]
      _end = buf[3]

      ithExonSize = _start - _prevEnd
      if (ithExonSize < 1){outflag = "ng"}
      blockSizes = blockSizes ( _start - _prevEnd ) ","

      # for (i+1)th exon
      blockStarts = blockStarts ( _end - chromStart ) ","
    }

    # --- exceptional cases ---
    if ( chromEnd < _end ) { chromEnd = _end + 1}
    # -------------------------

    blockSizes = blockSizes ( chromEnd - _end ) ","

    outLine = chrom "\t" chromStart "\t" chromEnd
    outLine = outLine "\t" name "\t" readsN "\t" strand
    outLine = outLine "\t" chromStart "\t" chromStart "\t" "0,0,0"
    outLine = outLine "\t" blockCount "\t" blockSizes "\t" blockStarts

    if (outflag == "ok") {
      print outLine
    } else {
      stde = tmpdir "/err.intron_readsTobed12.txt"
      printf "IncorrectModel: %s\n", outLine >> stde
    }
  }'
}


bed12ToBed12detail ()
{
  local prefix=$1
  awk --assign prefix=$prefix 'BEGIN{OFS="\t"}{
    buf = $4
    $4 = sprintf("%s%08d",prefix,NR)
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
  awk --assign minRefCount=$minRefCount --assign tmpdir=${tmpdir} 'BEGIN{OFS="\t"}{
    localNamesN=split($2,localNames,",")
    refNamesN=split($3,refNames,",")
    flag = "ok"
    for (i=1;i<=localNamesN;i++){
      str = sprintf("I%04d:%04d",i,localNamesN)
      if ( localNames[i] != str ) { flag = "ReadFilteredOut: " localNames[i] " != " str}
    }
    for (i=1;i<=refNamesN;i++){
      split( refNames[i] , buf, ":")
      refCount = buf[5]
      if ( refCount < minRefCount ) {flag = "ReadFilteredOut: count of " refNames[i] " is smaller than " minRefCount}
    }
    if (flag != "ok") {
      stde = tmpdir "/err.intronSetFilter.txt"
      printf "%s\t%s\n", flag, $0 >> stde
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
  | sort -k1,1 -k2,2n ${SORT_OPT_BASE} \
  > ${infile}.tmpboundary

  bedtools map -s -f $frac -r \
    -a ${infile}.tmpboundary -b ${infile}.tmpboundary \
    -c 2,3 -o mode,mode \
  | awk 'BEGIN{OFS="\t"}{
      boundaryMf = $1 ":" $7 ":" $8 ":" $6
      print $4, boundaryMf
    }'\
  | sort -k1,1 ${SORT_OPT_BASE} \
  > ${infile}.tmpboundaryMf

  join $infile ${infile}.tmpboundaryMf \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}'
}


bed12ToIntronCoordsName ()
{
  awk 'BEGIN{OFS="\t"}{
    chrom = $1; chromStart = $2; chromEnd = $3; 
    name = $4; score = $5; strand = $6; blockCount = $10;
    split($11,blockSizes,","); split($12,blockStarts,",");
    intronCoords = ""
    for (i = 2; i<= blockCount; i++)
    {
      intronStart = chromStart + blockStarts[i-1] + blockSizes[i-1]
      intronEnd = chromStart + blockStarts[i]
      intronCoord = sprintf("%s:%s:%s:%s", chrom, intronStart, intronEnd, strand)
      intronCoords = intronCoords intronCoord ","
    }
    if (intronCoords != "") {
      print intronCoords, name
    }
  }'
}



addIntronsMatchRef ()
{
  # annFile as bed12.gz
  local annFile=$1
  local tmpf=$( mktemp --tmpdir=${tmpdir} )

  awk '{print}' > ${tmpf}

  # prep intronCoordsName for input
  cat ${tmpf} \
  | bed12ToIntronCoordsName \
  | sort -k1,1 ${SORT_OPT_BASE} \
  > ${tmpf}_intronCoords_name.txt

  # prep intronCoordsName for annotation
  gunzip -c $annFile \
  | bed12ToIntronCoordsName \
  | sort -k1,1 ${SORT_OPT_BASE} \
  > ${tmpf}_intronCoords_nameAnn.txt

  # join
  join -t "	" -a 1 \
    ${tmpf}_intronCoords_name.txt \
    ${tmpf}_intronCoords_nameAnn.txt \
  | awk 'BEGIN{OFS="\t"}{if($3 == ""){$3 = "NA"}; print}' \
  | cut -f 2- \
  | sort -k1,1 ${SORT_OPT_BASE} \
  | groupBy -g 1 -c 2 -o collapse \
  | awk '{printf "%s\t%s,matchRef=%s\n", $1,$1,$2}' \
  > ${tmpf}_nameNewName.txt

  # print
  cat ${tmpf} \
  | awk 'BEGIN{OFS="\t"}{print $4, $0}' \
  | sort -k1,1 ${SORT_OPT_BASE} \
  | join -t "	" - ${tmpf}_nameNewName.txt \
  | cut -f 2- \
  | awk 'BEGIN{OFS="\t"}{$4 = $15;print}'  \
  | cut -f 1-14 \
  | sort -k1,1 -k2,2n ${SORT_OPT_BASE}
}



addIntronsMatchRefSingleExon ()
{
  # annFile as bed12.gz
  local annFile=$1
  local tmpf=$( mktemp --tmpdir=${tmpdir} )
  local tmpfa=$( mktemp --tmpdir=${tmpdir} )

  sort -k1,1 -k2,2n ${SORT_OPT_BASE} \
  > ${tmpf}

  gunzip -c $annFile \
  | awk '{if($10 == 1){print}}' \
  | sort -k1,1 -k2,2n ${SORT_OPT_BASE} \
  > ${tmpfa}

  # list ID pairs
  intersectBed -f 0.5 -r -nonamecheck -wa -wb -s -a ${tmpf} -b ${tmpfa} \
  | cut -f 4,18 \
  | sort -k1,1 -k2,2n ${SORT_OPT_BASE} \
  | groupBy -g 1 -c 2 -o collapse \
  | awk '{printf "%s\t%s,matchRef=%s\n", $1,$1,$2}' \
  > ${tmpf}_nameNewName.txt

  # add NA
  cut -f 4 ${tmpf} \
  | sort -k1,1 -k2,2n \
  | join -t "	" -a 1 - ${tmpf}_nameNewName.txt \
  | awk 'BEGIN{OFS="\t"}{if ($2 == ""){$2 = "NA"};print }' \
  > ${tmpf}_nameNewName.txt.tmp
  mv -f ${tmpf}_nameNewName.txt.tmp ${tmpf}_nameNewName.txt

  # include the matched information
  cat ${tmpf} \
  | awk 'BEGIN{OFS="\t"}{print $4, $0}' \
  | sort -k1,1 ${SORT_OPT_BASE} \
  | join -t "	" - ${tmpf}_nameNewName.txt \
  | cut -f 2- \
  | awk 'BEGIN{OFS="\t"}{$4 = $15;print}'  \
  | cut -f 1-14 \
  | sort -k1,1 -k2,2n ${SORT_OPT_BASE}
}



fivePrimeAttr ()
{
  awk \
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
    sigCountIncl = 0
    for (i=1;i<=readsN;i++)
    {
      match(reads[i],/\.5end-([ATGCN]*)-([ATGCN]*)/,buf)
      unmatch = buf[1]
      alignmentStart = buf[2]
      if ( isSig(unmatch) == "yes" ) {
        sigCount++
        sigCountIncl++
      } else {
        if ( match(reads[i],/\.5end--G/) ) { sigCountIncl++ }
      }
    }
    sigRatio = sigCount / readsN
    sigRatioIncl = sigCountIncl / readsN
    $4 = sprintf("%s,capSigCount=%d",$4,sigCount)
    $4 = sprintf("%s,capSigRatio=%.2f",$4,sigRatio)
    #$4 = sprintf("%s,capSigCountIncl=%d",$4,sigCountIncl)
    #$4 = sprintf("%s,capSigRatioIncl=%.2f",$4,sigRatioIncl)
    print

  }'
}


addLastExonOverlapWithOtherInternalExons ()
{
  local tmpf=$( mktemp --tmpdir=${tmpdir} )

  touch ${tmpf}.exon_internal
  touch ${tmpf}.exon_last
  awk '{print}' > ${tmpf}

  cat ${tmpf} \
  | awk 'BEGIN{OFS="\t"}{if($10 != 1){print}}' \
  | awk --assign tmpf=${tmpf} 'BEGIN{
      OFS="\t";
      outfile_internal=tmpf".exon_internal"
      outfile_last=tmpf".exon_last"
    }
    {
    chrom = $1; chromStart = $2; chromEnd = $3; 
    name = $4; score = $5; strand = $6; blockCount = $10;
    split($11,blockSizes,",")
    split($12,blockStarts,",")

    if ( strand == "+" ) {
      eStart = chromStart + blockStarts[blockCount]
      eEnd   = chromStart + blockStarts[blockCount] + blockSizes[blockCount]
      print chrom, chromStart, eStart, "internal", "0", strand > outfile_internal
      print chrom, eStart, eEnd, name, "0", strand > outfile_last
    } else {
      eStart = chromStart + blockStarts[1]
      eEnd   = chromStart + blockStarts[1] + blockSizes[1]
      print chrom, eStart, eEnd, name, "0", strand > outfile_last
      print chrom, eEnd, chromEnd, "internal", "0", strand > outfile_internal
    }
  }'

  sort -k1,1 -k2,2n ${SORT_OPT_BASE} ${tmpf}.exon_internal > ${tmpf}.exon_internal.tmp
  mv -f ${tmpf}.exon_internal.tmp ${tmpf}.exon_internal

  sort -k1,1 -k2,2n ${SORT_OPT_BASE} ${tmpf}.exon_last > ${tmpf}.exon_last.tmp
  mv -f ${tmpf}.exon_last.tmp ${tmpf}.exon_last

  intersectBed -nonamecheck -c -s -f 0.5 -a ${tmpf}.exon_last -b ${tmpf}.exon_internal \
  | cut -f 4,7 \
  > ${tmpf}.exon_internal_overlaps

  cat ${tmpf} \
  | awk 'BEGIN{OFS="\t"}{if($10 == 1){print}}' \
  | cut -f 1-6 | sort -k1,1 -k2,2n ${SORT_OPT_BASE} \
  | intersectBed -nonamecheck -c -s -f 0.5 -a - -b ${tmpf}.exon_internal \
  | cut -f 4,7 \
  >> ${tmpf}.exon_internal_overlaps

  sort -k1,1 ${SORT_OPT_BASE}  ${tmpf}.exon_internal_overlaps > ${tmpf}.exon_internal_overlaps.tmp
  mv -f ${tmpf}.exon_internal_overlaps.tmp ${tmpf}.exon_internal_overlaps

  cat ${tmpf} \
  | awk 'BEGIN{OFS="\t"}{print $4,$0}' \
  | sort -k1,1 ${SORT_OPT_BASE} \
  | join -t "	" - ${tmpf}.exon_internal_overlaps \
  | cut -f 2- \
  | awk 'BEGIN{OFS="\t"}{
    $4 = sprintf("%s,lastExonOverlapWithOtherInternalExons=%i", $4, $15)
    print
  }' \
  | cut -f 1-14
}



#m54284U_200720_151958/101714429/ccs.3endDown-AAAAAAAAAAAAAAATTAGC.3endPas--NOTFOUND.5end-CG-G:chr5:138002847:138019921:-
threePrimeAttr ()
{
  local minRatioA=$1

  cat ${tmpf} | awk --assign minRatioA=$minRatioA 'BEGIN{OFS="\t"}
  {
    lastExonOverlapWithOtherInternalExons=0
    match($4,/lastExonOverlapWithOtherInternalExons=([0-9]+)/,buf)
    if (RLENGTH >= 0) { lastExonOverlapWithOtherInternalExons = buf[1] }

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
    print
  }'
}



reAssignIds ()
{
  sort -k1,1 -k2,2n ${SORT_OPT_BASE}  \
  | awk --assign prefix=$prefix 'BEGIN{OFS="\t"}{
    match($4,"(SingleExon|MultiExon)[0-9]+,");
    if ( RLENGTH > 0 )
    {
      id = sprintf("%s%08d",prefix,NR)
      $4 = id "," substr($4,RLENGTH + 1)
      $13 = id
      print
    }
  }'
}



###
### main
###

mapQ=20
support_min_frac_intron=0.999
support_min_frac_boundary=0.95
annotation_reference=
threePrimeFilterMinRatioA=0.5
#threePrimeFilterMinInternalPrimingRatio=0.5
prefix=SG
debug_dir=

usage ()
{
  cat <<EOF

  $0 -i infile -g genome 
    [-a REFERENCE_TRANSCRIPT_MODEL.bed12.gz ]
    [-q mapQ(${mapQ})] 
    [-f support_min_frac_intron(${support_min_frac_intron})] 
    [-r support_min_frac_boundary(${support_min_frac_boundary})] 
    [-b threePrimeFilterMinRatioA(${threePrimeFilterMinRatioA})] 
    [-e threePrimeFilterMinInternalPrimingRatio(${threePrimeFilterMinInternalPrimingRatio})] 
    [-x prefix('${prefix}')] 
    [-d intermediate_file_directory_for_debug ] 

EOF
  exit 1
}

### handle options
while getopts i:g:a:q:f:r:c:d:b:e:x:d: opt
do
  case ${opt} in
  i) infile=${OPTARG};;
  g) genome=${OPTARG};;
  a) annotation_reference=${OPTARG};;
  q) mapQ=${OPTARG};;
  f) support_min_frac_intron=${OPTARG};;
  r) support_min_frac_boundary=${OPTARG};;
  b) threePrimeFilterMinRatioA=${OPTARG};;
  x) prefix=${OPTARG};;
  d) debug_dir=${OPTARG};;
  *) usage;;
  esac
done

if [ ! -n "${infile-}" ]; then usage; fi
if [ ! -n "${genome-}" ]; then usage; fi



if [ ! -n "${genome}.fai" ]; then
  samtools faidx ${genome}
fi

cat ${genome}.fai | cut -f 1,2 | sort -k1,1 ${SORT_OPT_BASE} > ${tmpdir}/genome.chrom_sizes
samtools view -bq ${mapQ} ${infile} > ${tmpdir}/infile.bam
touch ${tmpdir}/err.dummy.txt


bamAddDownstreamNucAsSeqNameSuffix ${tmpdir}/infile.bam ${genome} 20 \
| bamAddPolyASignalAsSeqNameSuffix \
| bamAdd5endAsSeqNameSuffix \
| bamToBed -bed12 -i stdin  \
| tee ${tmpdir}/infile.bed12 \
| bed12ToIntron \
| sort -k1,1 -k2,2n ${SORT_OPT_BASE} \
> ${tmpdir}/intron.bed

most_freq_intron ${tmpdir}/intron.bed $support_min_frac_intron \
| sort -k1,1 -k2,2n $SORT_OPT_BASE \
> ${tmpdir}/intron_mf.bed

intersectBed -nonamecheck -sorted -s -wa -wb \
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
| bed12ToBed12detail MultiExon \
> ${tmpdir}/outmodel.bedDetail

# single exons
cat ${tmpdir}/infile.bed12 \
| awk 'BEGIN{OFS="\t"}{if($10 == 1){print}}' \
| cut -f 1-6 \
| awk 'BEGIN{OFS="\t"}{$4=$4":"$1":"$2":"$3":"$6; print}' \
| sort -k1,1 -k2,2n  $SORT_OPT_BASE \
| mergeBed -s -c 4,5,6 -o collapse,count,distinct \
| awk 'BEGIN{OFS="\t"}{
    thickStart = $2; thickEnd = $2; itemRgb = "0,0,0";
    blockCount = 1; blockSizes = $3 - $2","; blockStarts = 0",";
    print $0, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
  }' \
| bed12ToBed12detail SingleExon \
> ${tmpdir}/outmodel.singleExon.bedDetail

if [ ! -n "${annotation_reference-}" ]; then

  # add single exons
  cat ${tmpdir}/outmodel.singleExon.bedDetail \
  >> ${tmpdir}/outmodel.bedDetail

  cat ${tmpdir}/outmodel.bedDetail \
  | sort -k1,1 -k2,2n ${SORT_OPT_BASE} \
  | fivePrimeAttr \
  | addLastExonOverlapWithOtherInternalExons \
  | threePrimeAttr $threePrimeFilterMinRatioA \
  | gzip -c --fast \
  > ${tmpdir}/outmodel_attr.bedDetail.gz

else

  # matching references
  cat ${tmpdir}/outmodel.bedDetail \
  | addIntronsMatchRef ${annotation_reference} \
  > ${tmpdir}/outmodel_ann.bedDetail

  # matching references (single exon)
  cat ${tmpdir}/outmodel.singleExon.bedDetail \
  | addIntronsMatchRefSingleExon ${annotation_reference} \
  >> ${tmpdir}/outmodel_ann.bedDetail

  cat ${tmpdir}/outmodel_ann.bedDetail
  | fivePrimeAttr \
  | addLastExonOverlapWithOtherInternalExons \
  | threePrimeAttr $threePrimeFilterMinRatioA \
  | gzip -c --fast \
  > ${tmpdir}/outmodel_attr.bedDetail.gz

fi

cat ${tmpdir}/err.*.txt >&2 

gunzip -c ${tmpdir}/outmodel_attr.bedDetail.gz \
| reAssignIds

if [ -n "${debug_dir-}" ]; then 
  mv ${tmpdir} ${debug_dir}
fi
