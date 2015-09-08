## 09.03.15
## GMB Wayne State University
## Makefile to align RNA-seq reads using STAR aligner
genomeindex=/wsu/home/groups/piquelab/data/RefGenome/STAR/
waspDir=/wsu/home/groups/piquelab/tools/WASP/
snpsDir=/wsu/home/groups/piquelab/gmb/SNPs/1KG/wasp_input/
fastqs=../fastqs
bamFolder=./

Qsub.ppn=12
Qsub.q=mmtxq
Qsub.N=bams

.SECONDARY:

%.Qsub:
	while ((`qme | wc -l` > 100)); do sleep 1; done
	touch $@
	echo "cd $${PWD}; make $*" | qsub -q $(Qsub.q) -l nodes=1:ppn=$(Qsub.ppn) -N $(Qsub.N) -o $@ -e $@.e
	sleep 1;	

## Workflow: align, quality filt, remove mapping bias, remove duplicates
all_bams: $(patsubst $(fastqs)/%_R1.fastq.gz, %.out.Qsub, $(wildcard $(fastqs)/*_R1.fastq.gz))
all_fq: $(patsubst %_quality.sort.bam, %_quality.remap.fq1.gz.Qsub, $(wildcard *_quality.sort.bam))
#all_clean: $(patsubst %_quality.remap.fq1.gz, %_clean.bam.Qsub, $(wildcard *_quality.remap.fq1.gz))

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
	samtools sort $*.out/$*.Aligned.out.bam $*.out/$*.Aligned.sort
	samtools index $*.out/$*.Aligned.sort.bam
	samtools idxstats $*.out/$*.Aligned.sort.bam > $*.out/$*.Aligned.sort.counts.txt

## Process and quality filter the aligned reads
%_merged.bam: 
	samtools merge $@ $*_L*.out/$*_L*.Aligned.sort.bam
	samtools index $@
	samtools view -c $@ > $*_merged_count.txt

%_quality.bam: %_merged.bam		
	samtools view -b1 -q10 $^ > $@
	samtools index $@
	samtools view -c $@ > $*_quality_count.txt
	ln -s $@ $*_quality.sort.bam

# %_clean.bam: %_quality.bam
# 	samtools view -b -q10 $^ | samtools.old rmdup - $@
# 	samtools index $@
# 	samtools view -c $@ > $*_clean_count.txt

## Use WASP to remove mapping bias from mapped allele-specific reads
%_quality.remap.fq1.gz: %_quality.bam
	python $(waspDir)/mapping/find_intersecting_snps.py -s -p $^ $(snpsDir)

#%_quality.remap.sorted.bam: %_quality.remap.fq1.gz 
#	bwa mem -t $(Qsub.ppn) $(genomeindex) $^ $*_quality.remap.fq2.gz | samtools view -Sb1 - > $*_quality.remap.bam
#	samtools sort $*_quality.remap.bam $*_quality.remap.sorted
#	samtools index $@

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
	samtools sort $*.remap.out/$*.Aligned.out.bam $*_quality.remap.sorted
	samtools index $@
	samtools idxstats $@ > $*_quality.remap.sorted.counts.txt

%_quality.remap.quality.bam: %_quality.remap.sorted.bam
	samtools view -b1 -q10 $^ > $@
	samtools index $@
	samtools view -c $@ > $*_quality.remap.quality.count.txt

%_quality.remap.keep.sorted.bam: %_quality.remap.quality.bam
	python $(waspDir)/mapping/filter_remapped_reads.py -p $*_quality.to.remap.bam $^ $*_quality.remap.keep.bam $*_quality.to.remap.num.gz
	samtools sort $*_quality.remap.keep.bam $*_quality.remap.keep.sorted
	samtools index $@

%_snpFilt.bam: %_quality.remap.keep.sorted.bam %_quality.keep.bam
	samtools merge $*_quality.remap.merged.bam $^
	samtools sort $*_quality.remap.merged.bam $*_snpFilt
	samtools index $@
	samtools view -c $@ > $*_snpFilt.count.txt

# %_clean.bam: %_snpFilt.bam
# 	python $(waspDir)/mapping//rmdup_pe.py $^ $*_rmdup.bam
# 	samtools sort $*_rmdup.bam $*_clean
# 	samtools index $@

## Temporary fix to adjust '1' vs 'chr1' issue. Long term: modify ref genome fasta file
%_clean.bam: %_clean.oldchr.bam
	samtools view $^ \
	| awk -v OFS='\t' '{ print $$1, $$2, "chr"$$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10, $$11, $$12, $$13, $$14, $$15 }' \
	| cat header.txt - \
	| samtools view - -b1 > $@
all_clean: $(patsubst %_clean.oldchr.bam, %_clean.bam, $(wildcard *_clean.oldchr.bam))

clean:
	rm *sorted*
	rm *merged*
	rm *quality*
	rm *clean*

#ll DswP9-HT10_quality.to.remap.bam DswP9-HT10_quality.remap.quality.bam DswP9-HT10_quality.remap.keep.bam DswP9-HT10_quality.to.remap.num.gz
