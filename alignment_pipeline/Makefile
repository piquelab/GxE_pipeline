## 09.03.15
## GMB Wayne State University
## Makefile to align RNA-seq reads using STAR aligner
## and allele-specific post-processing with WASP
## STAR: https://github.com/alexdobin/STAR
## WASP: https://github.com/bmvdgeijn/WASP

genomeindex=/wsu/home/groups/piquelab/data/RefGenome/STAR/
waspDir=/wsu/home/groups/piquelab/tools/WASP/
snpsDir=/wsu/home/groups/piquelab/gmb/SNPs/1KG/20130502/wasp_input/
fastqs=../fastqs/
bamFolder=./
Qsub.ppn=12
Qsub.q=mmtxq
Qsub.N=bams


## Utility commands & job scheduling
#.SECONDARY:

%.Qsub:
	while ((`qme | wc -l` > 100)); do sleep 1; done
	touch $@
	echo "cd $${PWD}; make $*" | qsub -q $(Qsub.q) -l nodes=1:ppn=$(Qsub.ppn) -N $(Qsub.N) -o $@ -e $@.e
	sleep 1;	

## Workflow: align, quality filt, remove mapping bias, remove duplicates
all_bams: $(patsubst $(fastqs)/%_R1.fastq.gz, %.out.Qsub, $(wildcard $(fastqs)/*_R1.fastq.gz))
all_quality: $(patsubst %_L1.out, %_quality.bam.Qsub, $(wildcard *_L1.out))
all_fq: $(patsubst %_quality.bam, %_quality.remap.fq1.gz.Qsub, $(wildcard *_quality.bam))
all_clean: $(patsubst %_quality.remap.fq1.gz, %_clean.bam.Qsub, $(wildcard *_quality.remap.fq1.gz))


## Align the sequencing reads using the STAR aligner
%.out: $(fastqs)/%_R1.fastq.gz $(fastqs)/%_R2.fastq.gz
	mkdir -p $@
	STAR \
		--runThreadN $(Qsub.ppn) \
		--genomeDir $(genomeindex) \
		--readFilesIn $^ \
		--readFilesCommand zcat \
		--outFileNamePrefix $*.out/$*. \
		--outSAMtype BAM Unsorted
	samtools sort -@ $(Qsub.ppn) $*.out/$*.Aligned.out.bam $*.out/$*.Aligned.sort
	samtools index $*.out/$*.Aligned.sort.bam
	samtools idxstats $*.out/$*.Aligned.sort.bam > $*.out/$*.Aligned_idxstats.txt

## Merge technical replicates and quality filter
%_merged.bam: 
	samtools merge -@ $(Qsub.ppn) $@ $*_L*.out/$*_L*.Aligned.sort.bam
	samtools index $@
	samtools view -c $@ > $*_merged_count.txt

%_quality.bam: %_merged.bam		
	samtools view -b1 -q10 $^ > $@
	samtools index $@
	samtools view -c $@ > $*_quality_count.txt
	ln -s $@ $*_quality.sort.bam

## Use WASP to identify potential mapping bias in reads overlapping SNPs
%_quality.remap.fq1.gz: %_quality.bam
	python $(waspDir)/mapping/find_intersecting_snps.py -s -p $^ $(snpsDir)

## Realign and quality filter reads overlapping SNPs 
%.remap.out: %_quality.remap.fq1.gz %_quality.remap.fq2.gz
	mkdir -p $@
	STAR \
		--runThreadN $(Qsub.ppn) \
		--genomeDir $(genomeindex) \
		--readFilesIn $^ \
		--readFilesCommand zcat \
		--outFileNamePrefix $*.remap.out/$*. \
		--outSAMtype BAM Unsorted

%_quality.remap.sorted.bam: %.remap.out
	samtools sort -@ $(Qsub.ppn) $*.remap.out/$*.Aligned.out.bam $*_quality.remap.sorted
	samtools index $@
	samtools idxstats $@ > $*_quality.remap.sorted_idxstats.txt

%_quality.remap.quality.bam: %_quality.remap.sorted.bam
	samtools view -b1 -q10 $^ > $@
	samtools index $@
	samtools view -c $@ > $*_quality.remap_count.txt

## Remove reads that map to different loci depending on the allele
%_quality.remap.keep.sorted.bam: %_quality.remap.quality.bam
	python $(waspDir)/mapping/filter_remapped_reads.py -p $*_quality.to.remap.bam $^ $*_quality.remap.keep.bam $*_quality.to.remap.num.gz
	samtools sort -@ $(Qsub.ppn) $*_quality.remap.keep.bam $*_quality.remap.keep.sorted
	samtools index $@

## Merge filtered reads and reads not overlapping SNPs
%_snpFilt.bam: %_quality.remap.keep.sorted.bam %_quality.keep.bam
	samtools merge -f -@ $(Qsub.ppn) $*_quality.remap.merged.bam $^ ## TODO: remove -f
	samtools sort -@ $(Qsub.ppn) $*_quality.remap.merged.bam $*_snpFilt
	samtools index $@
	samtools view -c $@ > $*_snpFilt_count.txt

## Remove duplicates
%_clean.bam: %_snpFilt.bam
	python $(waspDir)/mapping//rmdup_pe.py $^ $*_rmdup.bam
	samtools sort -@ $(Qsub.ppn) $*_rmdup.bam $*_clean
	samtools view -c $@ > $*_clean_count.txt
	samtools index $@

## TODO: add 'clean_intermediate' that will remove intermediate bam files
clean:
	rm *sorted*
	rm *merged*
	rm *quality*
	rm *clean*
