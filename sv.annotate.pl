#!/apps/perl/perl5/perls/perl-5.10.1/bin/perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use lib '/home/yanghao/perl5/lib/perl5/';
#use List::MoreUtils ':all';
use List::Util qw/uniq sum/;
use File::Basename;

my $base = dirname $0;

my ($vcf,$out) = @ARGV;

unless($vcf && -e $vcf && $out ){
	print my $eof = <<EOF;

=========================================================

speedseq sv 输出的sv结果（vcf文件）注释

用法：perl $0 <Illumina-B1701_sv.sv.vcf> <out.xls>

有问题请联系yanghao\@eulertechnology.com

=========================================================

EOF
	exit ;
}

my %sv_parser = (
	'DEL' => 'Deletion relative to the reference',
	'INS' => 'Insertion of novel sequence relative to the reference',
	'DUP' => 'Region of elevated copy number relative to the reference',
	'INV' => 'Inversion of reference sequence',
	'DUP' => 'Tandem duplication',
	'CNV' => 'Copy number variable region (may be both deletion and duplication)',
	'BND' => 'A breakend record is identified with the tag BND ref:https://samtools.github.io/hts-specs/VCFv4.2.pdf',
);
my %normal_chrom = (
			"chr1"=>undef,
			"chr2"=>undef,
			"chr3"=>undef,
			"chr4"=>undef,
			"chr5"=>undef,
			"chr6"=>undef,
			"chr7"=>undef,
			"chr8"=>undef,
			"chr9"=>undef,
			"chr10"=>undef,
			"chr11"=>undef,
			"chr12"=>undef,
			"chr13"=>undef,
			"chr14"=>undef,
			"chr15"=>undef,
			"chr16"=>undef,
			"chr17"=>undef,
			"chr18"=>undef,
			"chr19"=>undef,
			"chr20"=>undef,
			"chr21"=>undef,
			"chr22"=>undef,
			"chrX"=>undef,
			"chrY"=>undef,
);

my (%refseq,%chr_count);
open F,"/gpfs/users/yanghao/database/annotation/humandb//hg19_refGene.txt";
while(<F>){
	chomp;
	my @a = split /\t/;
	next unless $a[1] =~ /NM_/;
	next if $a[2] =~ /_/;
	$chr_count{$a[2]} ++;
	$refseq{$a[2]}{$chr_count{$a[2]}}{'s'} = $a[4];
	$refseq{$a[2]}{$chr_count{$a[2]}}{'e'} = $a[5];
	$refseq{$a[2]}{$chr_count{$a[2]}}{'g'} = $a[-4];
}
close F;
#print Dumper \%refseq;

my %hot_fusion;
open F,"$base/fusion.hotspots.txt";
while(<F>){
	chomp;
	my @a = split ;
	$hot_fusion{"$a[0]-$a[1]"} = undef;
	$hot_fusion{"$a[1]-$a[0]"} = undef;
}
close F;

open O,">$out" || dir $?;
say O "bio_marker\tleft_gene\tleft_chr\tleft_breakpoint\tright_gene\tright_chr\tright_breakpoint\tsv_type\tsv_desc\tsv_maf";
open F,"grep -v ^# $vcf|";
while(<F>){
	chomp;
	my @a = split /\t/;

	my ($sv_type) = $a[7] =~ m|SVTYPE=([^;]+)|;
	next if $sv_type eq 'INS';
	my ($sv_end) = $a[7] =~ m|END=([^;]+)|;

	my @k = split /:/,$a[-3];
	my @v = split /:/,$a[-2];
	my %info = ();
	@info{@k} = @v;

	my $left_chr = $a[0] =~ /chr/ ? $a[0] : "chr$a[0]";
	my $left_breakpoint = $a[1];
	my ($right_chr,$right_breakpoint);

	if($a[4] =~ /\<\S+\>/){
		$right_chr = $left_chr;
		$right_breakpoint = $sv_end;
	}else{
		($right_chr,$right_breakpoint) = $a[4] =~ /([^:\[\]]+):(\d+)/;
		$right_chr = $right_chr =~ /chr/ ? $right_chr : 'chr'.$right_chr;
	}

	next unless exists $normal_chrom{$left_chr} && exists $normal_chrom{$right_chr};
	my $sv_ratio = sprintf "%.3f%%",( (split /,/,$info{'AD'})[1] / sum(split /,/,$info{'AD'}) )* 100;

	#左基因注释
	my $left_gene;
	for my $i (sort {$refseq{$left_chr}{$a} <=> $refseq{$left_chr}{$b}} keys %{$refseq{$left_chr}}){
		$left_gene .= "$refseq{$left_chr}{$i}{'g'};" if $left_breakpoint >= $refseq{$left_chr}{$i}{'s'} && $left_breakpoint <= $refseq{$left_chr}{$i}{'e'};
	}
	$left_gene ||= '-';
	$left_gene = join ';',uniq(split /;/,$left_gene);
	

	#右基因注释
	my $right_gene;
	for my $i (sort {$refseq{$right_chr}{$a} <=> $refseq{$right_chr}{$b}} keys %{$refseq{$right_chr}}){
		$right_gene .= "$refseq{$right_chr}{$i}{'g'};" if $right_breakpoint >= $refseq{$right_chr}{$i}{'s'} && $right_breakpoint <= $refseq{$right_chr}{$i}{'e'};
	}
	$right_gene ||= '-';
	$right_gene = join ';',uniq(split /;/,$right_gene);
	
	my $bio_marker = exists $hot_fusion{"$left_gene-$right_gene"}  ? "$left_gene-$right_gene" : exists $hot_fusion{"$right_gene-$left_gene"} ? "$right_gene-$left_gene" : '-';
	
	say O join "\t",($bio_marker,$left_gene,$left_chr,$left_breakpoint,$right_gene,$right_chr,$right_breakpoint,$sv_type,$sv_parser{$sv_type},$sv_ratio);
}
close F;
