use strict;
use warnings;
use lib '/home/yanghao/perl5/lib/perl5/';
use Data::Dumper;
use feature 'say';

my ($sample,$bam,$interval) = @ARGV;


my $eof = <<EOF;

perl $0 <sample> <bam> <interval.bed>

统计bam文件的各种参数，例如覆盖度、平均测序深度等指标

EOF


die $eof unless $sample && $bam && -e $bam && $interval && -e $interval;


chomp(my $total_reads = `head -n1 $bam.flagstat|cut -d' ' -f1`);
chomp(my $on_target_reads = `sambamba view -t\`nproc\` -L $interval $bam -c`);
my $cap_ratio = sprintf "%.3f",$on_target_reads/$total_reads;
my $total_base = $total_reads * 150;
my $total_base_m = $total_base / 1000000;
my $total_base_g = $total_base_m / 1000;
chomp(my $dup_reads = `grep duplicates $bam.flagstat |cut -d' '  -f1`);
my $dup_ratio = $dup_reads / $total_reads;
chomp(my $map_ratio = `grep mapped $bam.flagstat|head -n1|perl -lne'\$_=~m|\\((.*\\%):|;print \$1'`);
chomp(my $average_depth = `zcat $sample.stat.regions.bed.gz|perl -lane'\$s+=\$F[-1]}{print \$s/\$.'`);
my ($target_size,$cov_1,$cov_10,$cov_20,$cov_30,$cov_100,$cov_300);
open F,"zcat $sample.stat.thresholds.bed.gz|sed 1d|";
while(<F>){
	chomp;
	my @a = split;
	$target_size += $a[2] - $a[1];
	$cov_1 += $a[4];
	$cov_10 += $a[5];
	$cov_20 += $a[6];
	$cov_30 += $a[7];
	$cov_100 += $a[8];
	$cov_300 += $a[9];
}
close F;
$cov_1 /= $target_size;
$cov_10 /= $target_size;
$cov_20 /= $target_size;
$cov_30 /= $target_size;
$cov_100 /= $target_size;
$cov_300 /= $target_size;

print "sample\tbam\ttotal_reads\ttotal_base(bp)\ttotal_base(M)\ttotal_base(G)\tdup_reads\tdup_ratio\tmap_ratio\taverage_depth_on_target\ttarget_size\tgt_1x_on_target\tgt_10x_on_target\tgt_20x_on_target\tgt_30x_on_target\tgt_100x_on_target\tgt_300x_on_target\ton_target_reads\ton_target_reads_ratio\n";
print "$sample\t$bam\t$total_reads\t$total_base\t$total_base_m\t$total_base_g\t$dup_reads\t$dup_ratio\t$map_ratio\t$average_depth\t$target_size\t$cov_1\t$cov_10\t$cov_20\t$cov_30\t$cov_100\t$cov_300\t$on_target_reads\t$cap_ratio\n";
