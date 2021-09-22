
# spectra - building gene models based on full-length cDNA long reads

spectra (SPlice matching, End Consistent TRAnscripts) is a set of scripts
to build gene models based on full-length cDNA reads. It employs an approach
to enrich 5'-complete reads by using 5'-end nucleotides not encoded in the genome.

It has to be noted that preference of unencoded nucleotides at 5-ends largely depends
on a protocol to produce cDNAs. This alpha version is *tentatively* optimized for the protocol
based on template-switch [^1], sequenced by PacBio CCS. `pbclean.sh`, a helper script, is
a simple implementation of ENCODE Long Read RNA-Seq Analysis Protocol for Human Samples (v3.1) [^2],
with some modification to handle 5'-end carefully. 


Requirements
-------------
* [minimap2](https://github.com/lh3/minimap2) tested in 2.22
* [bedtools](https://github.com/arq5x/bedtools2) tested in 2.29.2
* [samtools](https://github.com/samtools/samtools) tested in v1.13
* [jksrc](https://hgdownload.soe.ucsc.edu/admin/) tested in v385
* [flip_reads.py](https://github.com/mortazavilab/ENCODE-references) tested in f3b387e commit 
* [isoseq3](https://github.com/PacificBiosciences/IsoSeq) tested in v3.4.0
* [lima](https://github.com/PacificBiosciences/barcoding) tested in v2.2.0


Steps to run
------------
1. align fastq by `twoPassAlign.sh` (input/output format: FASTQ/BAM)
    - it finds relatively frequent junctions in the first pass, which is prioritized in the second pass
    - two pass alignments is employed to avoid small inconsistencies caused by sequencing error
    - when reference gene models ( genemodel.bed12.gz ) is specified, their junctions are prioritized in the first pass
2. group alignments into isoform group by `speGroup.sh` (input/output format: BAM/BED12+2[that is, bedDetail])
    - alignments sharing the same set of introns are grouped
    - features on their 5'-end / 3'-end alignments are added to both individual reads and isoform groups
3. select representative models by `speSelect.sh` (input/output format: BED12+2, BED12+2)


Helper scripts
--------------
* `speFilter.sh` - used in `speSelect.sh`
* `gtf2bed12.sh` - format conversion, which is used to convert reference gene model into BED12
* `pbclean.sh` - take CCS raw reads in BAM to produce clean CCS reads in FASTQ. 5'-ends are carefully processed.
* `seqnamesuffix2fq.sh` - add suffix to sequence names, which can be used to add replicate information


End features (inspected by `speGroup.sh`)
------------------------------------------
* 5-end cap signature
    - `5end-{G|C}*G[GENOME_ALIGNMENT]`, where the number of `G` should be larger than `C`
* Internal priming signature - when both of the following conditions are true
    - No major PAS (poly adenylation signal,  `AAUAAA` or `AUUAAA`)
    - Downstream (20bp) of the alignments on the genome is `A` rich (>50%)


Selection criteria (by `speSelect.sh`)
---------------------------------------
* matching to the reference 
* 5-end cap signature (counts, ratio) - models with this signatures are selected
    - models with only one read having this signature is less prioritize, in case that other overlapping models have more
* internal priming - models meeting both of the following conditions are filtered out.
    - Internal priming signature (counts, ratio)
    - the last exon overlap other gene models of their non last exon
* at least one gene model, per locus.


Citation
------

    Kawaji H. (2021) spectra, a set of scripts to build gene models based on full-length cDNA reads, https://github.com/hkawaji/spectra

---

[^1]: Protocol to build non-size selected cDNA libraries for Pacific Biosciences
    long-read sequencing Version 3.0 (October, 2020)
    https://www.encodeproject.org/experiments/ENCSR507JOF/

[^2]: ENCODE Long Read RNA-Seq Analysis Protocol for Human Samples (v3.1)
  https://www.encodeproject.org/experiments/ENCSR634AKC/
