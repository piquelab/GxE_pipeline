## 3.13.14
## rpique [2014-03-13  8:07 ] I changed the resource requirments for sorting. 
## gmoyerbr [2015-08-21 1:56 ] Integrated WASP-based postprocessing
## Make .bams for system plates
genomeindex=/wsu/home/groups/piquelab/data/RefGenome/hg19.fa.gz
fastqs=../fastqs
waspDir=/wsu/home/groups/piquelab/tools/WASP/
sortedBams=$(shell find ./ -name '*sorted.bam' | sed 's/_[SL].*\.bam//g' | sort | uniq)
Qsub.ppn=2
Qsub.q=mmtxq
Qsub.N=bams

## Utility commands & job scheduling
all: 
	echo $(sortedBams)
%.Qsub:
	while ((`qme | wc -l` > 100)); do sleep 1; done
	touch $@
	echo "cd $${PWD}; make $*" | qsub -q $(Qsub.q) -l nodes=1:ppn=$(Qsub.ppn) -N $(Qsub.N) -o $@ -e $@.e
	sleep 1;	


## Initial alignment of FASTQ files to specified genome
%.bam: $(fastqs)/%_R1.fastq.gz $(fastqs)/%_R2.fastq.gz
	bwa mem -t $(Qsub.ppn) $(genomeindex) $^  | samtools view -Sb1 - > $@
	samtools sort $@ $*_sorted
	samtools index $*_sorted.bam

all_bams: $(patsubst $(fastqs)/%_R1.fastq.gz,%.bam.Qsub,$(wildcard $(fastqs)/*_R1.fastq.gz))


## Merge technical replicates and quality filter
%_merged.bam: %_[SL]*_sorted.bam 
	samtools merge $@ $^
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
# all_clean: $(patsubst %, %_clean.bam.Qsub,$(sortedBams))


## Use WASP to remove mapping bias from mapped allele-specific reads
%_quality.remap.fq1.gz: %_quality.bam
	python $(waspDir)/mapping/find_intersecting_snps.py -p $^ SNP_file_directory

%_sorted_mapw.bam: %_quality.remap.fq1.gz 
	bwa mem -t $(Qsub.ppn) $(genomeindex) $^ $*_quality.remap.fq2.gz | samtools view -Sb1 - > $*_mapw.bam
	samtools sort $*_mapw.bam $*_sorted_mapw
	samtools index $@

%_quality_mapw.bam: %_sorted_mapw.bam
	samtools view -b1 -q10 $^ > $@
	samtools index $@
	samtools view -c $@ > $*_quality_mapw_count.txt

%_quality.remap.keep.sorted.bam: %_quality_mapw.bam 
	python $(waspDir)/mapping//filter_remapped_reads.py -p $*_quality.to.remap.bam $^ $*_quality.remap.keep.bam $*_quality.to.remap.num.gz
	samtools sort $*_quality.remap.keep.bam $*_quality.remap.keep.sorted
	samtools index $@

%_.clean.bam: %_quality.remap.keep.sorted.bam
	python $(waspDir)/mapping//rmdup_test.py $^ $@

clean:
	rm *sorted*
	rm *merged*
	rm *quality*
	rm *mapw*
	rm *remap*
	rm *clean*
