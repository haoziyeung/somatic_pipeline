#!/apps/perl/perl5/perls/perl-5.10.1/bin/perl
use strict;
use warnings;

use lib '/home/yanghao/perl5/lib/perl5/';
use Term::ANSIColor;
use Getopt::Long qw/GetOptions/;
use Data::Dumper;
use Cwd;
use File::Basename;

my $base = dirname $0;
use feature 'say';

my ($tumor,$normal,$t1,$t2,$n1,$n2,$bed,$part);
my $eof = <<EOF;

----------------------------------------------------------
	{体细胞突变检测流程}
特点：
1)使用sentieon TNscope流程，准确、高效、快速
2)关联cosmic、www.cancergenomeinterpreter.org靶药数据库
3)分析SV，常见的融合基因一目了然
4)分析CNV，可检测拷贝数变异


输入：fastq

输出：somatic_mut.xlsx/tsv,SV,CNV,QC metrics...

使用方法：
	perl $0 -t <tumor name> -t1 <tumor R1> -t2 <tumor R2> -n <normal name> -n1 <normal R1> -n2 <normal R2> -b <interval.bed>

默认提交到cn-medium分区，如果没有资源了，你可以指定作业提交到cn-fast，如下所示：
	perl $0 -t <tumor name> -t1 <tumor R1> -t2 <tumor R2> -n <normal name> -n1 <normal R1> -n2 <normal R2> -b <interval.bed> -part cn-fast

Last Updated:2017/12/25

如果使用过程中发生任何问题，请联系yanghao\@eulertechnology.com

---------------------------------------------------------

EOF

GetOptions(
    't=s' => \$tumor,
    't1=s' => \$t1,
    't2=s' => \$t2,

    'n=s' => \$normal,
    'n1=s' => \$n1,
    'n2=s' => \$n2,

    'b=s' => \$bed,
    'part=s' => \$part,
);

unless($tumor && $t1 && -e $t1 && $t2 && -e $t2 && $normal && $n1 && -e $n1 && $n2 && -e $n2){
    say colored($eof, "bright_yellow");
    exit;
}


$bed ||= '/gpfs/users/yanghao/project/shi-jian-zhi-ping/t_n_20171030/illumina/target.bed';
$part ||= 'cn-medium';
my @partition = qw/cn-fast cn-medium/;
die "只能提交到SLURM作业调度系统的cn-fast或者cn-medium分区上\n" unless $part ~~ @partition;

my $account = $part eq 'cn-medium' ? 'cnm' : 'cnf';
my $path = getcwd;
my $ref = '/gpfs/users/yanghao/database/ref/b37/human_g1k_v37_decoy.fasta';
my $sentieon = 'export SENTIEON_LICENSE=io03:8990;/gpfs/bin/sentieon/bin/sentieon';
my $nt = '`nproc`';
my $known_Mills_indels = '/gpfs/users/yanghao/database/gatk_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz';
my $known_1000G_indels = '/gpfs/users/yanghao/database/gatk_bundle/1000G_phase1.indels.b37.vcf.gz';
my $dbsnp = '/gpfs/users/yanghao/database/ncbi/dbsnp/All_20170710.vcf.gz';

# ******************************************************************************
# TUMOR fastq alignment,indel realignment,base recalibration
# ******************************************************************************
open F,">$path/somatic.$tumor.align.sh";
print F <<EOF;
#!/bin/bash

#SBATCH -p $part
#SBATCH -N 1
#SBATCH --no-requeue
#SBATCH -A $account
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15

cd $path

. /home/yanghao/.bashrc

speedseq align -t $nt -o $tumor -R \"\@RG:\\tID:$tumor\\tSM:$tumor\\tPL:ILLUMINA\\tCN:Euler\\tLB:test_lb\\tPU:test_pu\" $ref $t1 $t2

sambamba flagstat -t $nt $tumor.bam > $tumor.bam.flagstat

/home/yanghao/anaconda2/pkgs/mosdepth-0.2.0-htslib1.6_0/bin/mosdepth -t 4 -b $bed -n -T1,10,20,30,100,300 $tumor.stat $tumor.bam

$sentieon driver -r $ref -t $nt -i $tumor.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels $tumor.realigned.bam

$sentieon driver -r $ref -t $nt -i $tumor.realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels $tumor.recal_data.table

EOF
close F;
chomp(my $tumor_aln= `sbatch $path/somatic.$tumor.align.sh`);
$tumor_aln =~ s/Submitted batch job //g;

# ******************************************************************************
# NORMAL fastq alignment,indel realignment,base recalibration
# ******************************************************************************
open F,">$path/somatic.$normal.align.sh";
print F <<EOF;
#!/bin/bash

#SBATCH -p $part
#SBATCH -N 1
#SBATCH --no-requeue
#SBATCH -A $account
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15

cd $path

. /home/yanghao/.bashrc

speedseq align -t \`nproc\` -o $normal -R \"\@RG:\\tID:$normal\\tSM:$normal\\tPL:ILLUMINA\\tCN:Euler\\tLB:test_lb\\tPU:test_pu\" $ref $n1 $n2

sambamba flagstat -t\`nproc\` $normal.bam > $normal.bam.flagstat

/home/yanghao/anaconda2/pkgs/mosdepth-0.2.0-htslib1.6_0/bin/mosdepth -t 4 -b $bed -n -T1,10,20,30,100,300 $normal.stat $normal.bam

$sentieon driver -r $ref -t $nt -i $normal.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels $normal.realigned.bam

$sentieon driver -r $ref -t $nt -i $normal.realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels $normal.recal_data.table

EOF
close F;
chomp(my $normal_aln= `sbatch $path/somatic.$normal.align.sh`);
$normal_aln =~ s/Submitted batch job //g;

# ****************************************************************************
# Corealignment of tumor and normal,Somatic and Structural variant calling
# ****************************************************************************
open F,">$path/somatic.$tumor-$normal.corealn.sh";
print F <<EOF;
#!/bin/bash

#SBATCH -p $part
#SBATCH -N 1
#SBATCH --no-requeue
#SBATCH -A $account
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --dependency=$tumor_aln,$normal_aln

cd $path

$sentieon driver -r $ref -t $nt -i $tumor.realigned.bam -i $normal.realigned.bam -q $tumor.recal_data.table -q $normal.recal_data.table --algo Realigner -k $known_Mills_indels -k $known_1000G_indels $tumor-$normal.tn_corealigned.bam

$sentieon driver -r $ref -t $nt -i $tumor-$normal.tn_corealigned.bam --interval $bed --interval_padding 100 --algo TNscope --tumor_sample $tumor --normal_sample $normal --dbsnp $dbsnp --pcr_indel_model conservative --sv_mask_ext 10 --max_fisher_pv_active 0.05 --min_tumor_allele_frac 0.0005 --no_mapq_cap 1 --clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 1.0 --assemble_mode 4 --cosmic /gpfs/users/yanghao/database/cosmic/CosmicCodingMuts.vcf.gz --cosmic /gpfs/users/yanghao/database/cosmic/CosmicNonCodingVariants.vcf.gz $tumor-$normal.vcf.gz

zgrep -v ^# $tumor-$normal.vcf.gz | grep SVTYPE > $tumor-$normal.sv.vcf

zcat $tumor-$normal.vcf.gz | perl -lne'if(/^#/){print}else{print unless /SVTYPE/}' > $tumor-$normal.snp_indel.vcf

/gpfs/users/yanghao/software/annovar/table_annovar.pl $tumor-$normal.snp_indel.vcf /gpfs/users/yanghao/database/annotation/humandb/ -buildver hg19 -protocol refGene,ensGene,cytoBand,esp6500siv2_all,1000g2015aug_all,avsnp150,gnomad_genome,gnomad_exome,dbnsfp33a,cosmic70 -operation g,g,r,f,f,f,f,f,f,f -argument \"-splicing_threshold 5 --hgvs,,,,,,,,,\" -vcfinput -nastring . --thread 6 --maxgenethread 6 -polish

perl $base/parse_hgvs.pl $tumor-$normal.snp_indel.vcf.hg19_multianno.txt > $tumor-$normal.snp_indel.xls

perl $base/sv.annotate.pl $tumor-$normal.sv.vcf $tumor-$normal.sv.txt

perl $base/bam.qc.pl $tumor $tumor.bam $bed > $tumor.stat.txt

perl $base/bam.qc.pl $normal $normal.bam $bed > $normal.stat.txt

samtools bedcov $bed $tumor.bam > $tumor.bam.bedcov

Rscript /gpfs/bin/EulerCNV/MSKCC.b37.plot.CNVgraph.R $tumor.bam.bedcov

EOF
close F;
chomp(my $co_aln = `sbatch $path/somatic.$tumor-$normal.corealn.sh`);
$co_aln =~ s/Submitted batch job //g;

say "========================================\n作业提交成功\n=======================================\n";
print my $slurm = `squeue -lu$ENV{'USER'}`;
