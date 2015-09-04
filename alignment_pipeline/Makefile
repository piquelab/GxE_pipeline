## 09.03.15
## GMB Wayne State University
## Makefile to align RNA-seq reads using STAR aligner
genomeindex=/wsu/home/groups/piquelab/data/RefGenome/STAR/
fastqs=../fastqs
bamFolder=./

Qsub.ppn=2
Qsub.q=mmtxq
Qsub.N=bams


all: 
	echo $(sortedBams)

%.Qsub:
	while ((`qme | wc -l` > 100)); do sleep 1; done
	touch $@
	echo "cd $${PWD}; make $*" | qsub -q $(Qsub.q) -l nodes=1:ppn=$(Qsub.ppn) -N $(Qsub.N) -o $@ -e $@.e
	sleep 1;	

all_bams: $(patsubst $(fastqs)/%_R1.fastq.gz,%.bam.Qsub,$(wildcard $(fastqs)/*_R1.fastq.gz))

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

# %_clean.bam: %_quality.bam
# 	samtools view -b -q10 $^ | samtools.old rmdup - $@
# 	samtools index $@
# 	samtools view -c $@ > $*_clean_count.txt

## Use WASP to remove mapping bias from mapped allele-specific reads
%_quality.remap.fq1.gz: %_quality.bam
	ln -s $*_quality.bam $*_quality.sort.bam
	python $(waspDir)/mapping/find_intersecting_snps.py -s -p $^ $(snpsDir)

%_quality.remap.sorted.bam: %_quality.remap.fq1.gz 
	bwa mem -t $(Qsub.ppn) $(genomeindex) $^ $*_quality.remap.fq2.gz | samtools view -Sb1 - > $*_quality.remap.bam
	samtools sort $*_quality.remap.bam $*_quality.remap.sorted
	samtools index $@

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

%_clean.bam: %_snpFilt.bam
	python $(waspDir)/mapping//rmdup_pe.py $^ $*_rmdup.bam
	samtools sort $*_rmdup.bam $*_clean
	samtools index $@

clean:
	rm *sorted*
	rm *merged*
	rm *quality*
	rm *clean*
