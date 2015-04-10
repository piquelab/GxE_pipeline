sortedBams=$(shell find ./ -name '*quality.bam' | sed 's/_quality\.bam//g' )
genomeindex=/wsu/home/groups/piquelab/data/RefGenome/hg19.fa.gz
source=~/piquelab/charvey/WASP/mapping
Qsub.ppn=4
Qsub.q=mmtxq
Qsub.N=keepbams

.SECONDARY:

all_keep: $(patsubst %, %_quality.remap.keep.clean.bam.Qsub, $(sortedBams))

%.Qsub:
	touch $@
	echo "cd $${PWD}; module load python-2.7; make $*" | qsub -q $(Qsub.q) -l nodes=1:ppn=$(Qsub.ppn) -N $(Qsub.N) -o $@ -e $@.e
	sleep 1;

%_quality.remap.fq1.gz: %_quality.bam
	python $(source)/find_intersecting_snps.py -p $^ SNP_file_directory

%_sorted_mapw.bam: %_quality.remap.fq1.gz 
	bwa mem -t $(Qsub.ppn) $(genomeindex) $^ $*_quality.remap.fq2.gz | samtools view -Sb1 - > $*_mapw.bam
	samtools sort $*_mapw.bam $*_sorted_mapw
	samtools index $@

%_quality_mapw.bam: %_sorted_mapw.bam
	samtools view -b1 -q10 $^ > $@
	samtools index $@
	samtools view -c $@ > $*_quality_mapw_count.txt

%_quality.remap.keep.sorted.bam: %_quality_mapw.bam 
	python $(source)/filter_remapped_reads.py -p $*_quality.to.remap.bam $^ $*_quality.remap.keep.bam $*_quality.to.remap.num.gz
	samtools sort $*_quality.remap.keep.bam $*_quality.remap.keep.sorted
	samtools index $@

%_quality.remap.keep.clean.bam: %_quality.remap.keep.sorted.bam
	python $(source)/rmdup_test.py $^ $@
