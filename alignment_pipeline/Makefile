## 3.13.14
## rpique [2014-03-13  8:07 ] I changed the resource requirments for sorting. 
## gmoyerbr [2015-08-21 1:56 ] Integrated WASP-based postprocessing
## Make .bams for system plates
genomeindex=/wsu/home/groups/piquelab/data/RefGenome/hg19.fa.gz
fastqs=../fastqs
waspDir=/wsu/home/groups/piquelab/tools/WASP/
snpsDir=/wsu/home/groups/piquelab/gmb/SNPs/1KG/wasp_input/
sortedBams=$(shell find ./ -name '*sorted.bam' | sed 's/_[SL].*\.bam//g' | sort | uniq)
Qsub.ppn=2
Qsub.q=mmtxq
Qsub.N=bams


## Utility commands & job scheduling
.PRECIOUS: %_quality.bam

all: 
	echo $(sortedBams)

%.Qsub:
	while ((`qstat $(Qsub.q) | grep ' Q ' | wc -l` > 100)); do sleep 1; done
	touch $@
	echo "cd $${PWD}; make $*" | qsub -q $(Qsub.q) -l nodes=1:ppn=$(Qsub.ppn) -N $(Qsub.N) -o $@ -e $@.e
	sleep .5


## Workflow: align, quality filt, remove mapping bias, remove duplicates
all_bams: $(patsubst $(fastqs)/%_R1.fastq.gz, %.bam.Qsub, $(wildcard $(fastqs)/*_R1.fastq.gz))
all_quality: $(patsubst %, %_quality.bam.Qsub, $(sortedBams))
all_fq: $(patsubst %_quality.bam, %_quality.remap.fq1.gz.Qsub, $(wildcard *_quality.bam))
all_clean: $(patsubst %_quality.remap.fq1.gz, %_clean.bam.Qsub, $(wildcard *_quality.remap.fq1.gz))


## Initial alignment of FASTQ files to specified genome
## Merge technical replicates and quality filter
%.bam: $(fastqs)/%_R1.fastq.gz $(fastqs)/%_R2.fastq.gz
	bwa mem -t $(Qsub.ppn) $(genomeindex) $^  | samtools view -Sb1 - > $@
	samtools sort $@ $*_sorted
	samtools index $*_sorted.bam

%_merged.bam: %_[SL]*_sorted.bam 
	samtools merge $@ $^
	samtools index $@
	samtools view -c $@ > $*_merged_count.txt

%_quality.bam: %_merged.bam		
	samtools view -b1 -q10 $^ > $@
	samtools index $@
	samtools view -c $@ > $*_quality_count.txt


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


## Bleach the working directory
clean:
	rm *sorted*
	rm *merged*
	rm *quality*
	rm *remap*
	rm *snpFilt*
	rm *clean*
