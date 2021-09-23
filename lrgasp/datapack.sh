#!/bin/bash

#infile=out/mouse-ES-cDNA-PacBio-ENCSR172GXL.all.bedDetail.gz
#outprefix=abc

infile=$1
outprefix=$2

usage ()
{
  cat <<EOF
usage: $0 INFILE.bedDetail.gz OUTPREFIX

EOF
  exit 1
}

if [ ! -n "${infile-}" ]; then usage; fi
if [ ! -n "${outprefix-}" ]; then usage; fi



outfile_gtf_gz=${outprefix}.gtf.gz
outfile_read_model_map_gz=${outprefix}.tsv.gz

gunzip -c out/mouse-ES-cDNA-PacBio-ENCSR172GXL.all.bedDetail.gz \
| cut -f 13,14 \
| awk 'BEGIN{OFS="\t"}{
    transcript_id = $1
    n = split($2,reads,",")
    for (i = 1; i <= n; i++)
    {
      print $1, reads[i], transcript_id
    }
  }' \
| cut -f 1-3 -d '/' \
| awk 'BEGIN{
    OFS="\t"
    printf "read_id\ttranscript_id\n"
  }{
    print $2, $1
  }' \
| gzip -c \
> ${outfile_read_model_map_gz}


gunzip -c out/mouse-ES-cDNA-PacBio-ENCSR172GXL.all.bedDetail.gz \
| awk 'BEGIN{OFS="\t"}{
    id = "."; matchRef= "NA"
    match($4, "^[^,]+")
    if ( RLENGTH > 0 ) { id = substr($4, RSTART, RLENGTH) }
    match($4, "matchRef=[^,]+")
    if ( RLENGTH > 0 ) { matchRef = substr($4, RSTART + 9, RLENGTH - 9) }
    $4 = id "," matchRef
    #$4 = id
    print
  }' \
| cut -f 1-12 \
| bedToGenePred stdin stdout \
| genePredToGtf file stdin stdout \
| awk 'BEGIN{FS="\t";OFS="\t"}{
    if ( $3 != "exon") {next}

    #exon_number
    match($9, "exon_number \"[0-9]+\";")
    en = substr($9, RSTART, RLENGTH)

    #ids
    match($9, "\"[^,]+,[^\"]+\"")
    if ( RLENGTH > 0 ) { buf = substr($9, RSTART , RLENGTH ) }
    gsub("\"","",buf); split(buf,ids,",");
    id = ids[1]; refid = ids[2];
    name = sprintf("gene_id \"%s\"; transcript_id \"%s\"; %s", id, id, en)
    if ( refid != "NA" ) { name = name " reference_transcript_id \"" refid "\";"}
    $9 = name
    print
  }' \
| gzip -c \
> ${outfile_gtf_gz}


