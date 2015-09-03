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
%.merged.bam: 
	samtools merge $@ $*_L*.out/$*_L*.Aligned.sort.bam
	samtools index $@
	samtools view -c $@ > $*_merged_count.txt

%_quality.bam: %_merged.bam		
	samtools view -b1 -q10 $^ > $@
	samtools index $@
	samtools view -c $@ > $*_quality_count.txt

%_clean.bam: %_quality.bam
	samtools view -b -q10 $^ | samtools.old rmdup - $@
	samtools index $@
	samtools view -c $@ > $*_clean_count.txt
clean:
	rm *sorted*
	rm *merged*
	rm *quality*
	rm *clean*
