#!usr/bin/perl

#Copyright  (c) 2022  Soichirou Satoh
#This Perl script is released under the MIT License.
#http://opensource.org/licenses/mit-license.php

# Automatically filtering of all paired-end fastq files in the current directory using Trimmomatic.

##### I. Set arguments ###################################################################
use Cwd 'getcwd';
$curdir = getcwd;
$outdir = $curdir . "\/Fastq_output";

use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);
$adapterfa = ""; # -a
$minlen = 50; # -m
$threadnum = 1; # -t
$help = ""; # -h
GetOptions(
	'a=s' => \$adapterfa,
	'm=i' => \$minlen,
	't=i' => \$threadnum,
	'h' => \$help,
) || die "Invalid options";

unless ($help == 1) {
	print "Automatically filtering of all paired\-end fastq files \(_1\.fastq\, _2\.fastq\) in the current directory\.\n";
	print "Proceed \(y\/n\)\?\: ";
	$check = <STDIN>;
	chomp $check;
	$check = lc($check);
	unless ($check eq "y") {
		$help = 1;
	}
}

if ($help == 1) {
	print "\nFilteringFastq_for_MPASS\.pl\n\n";
	
	print "Usage\: perl FilteringFastq_for_MPASS\.pl \[options\]\n\n";
	print "Options\:\n";
	print "\t\-a Adapter fasta file \[String\]\n";
	print "\t\-m Minimum length of filtered fastq reads \(Default\: 50\) \[Integer\]\n";
	print "\t\-t Number of threads \(Default\: 1\) \[Integer\]\n\n";
	print "\t\-h Help\n\n\n";

	print "Requirements\:\n";
	print "\tPaired-end fastq files \(_1\.fastq\, _2\.fastq\)\n\n";
	print "\t\[Software \(Commands and URL for installation\)\]\n";
	print "\tTrimmomatic \(conda install \-c bioconda trimmomatic\)\n";
	print "\tAdapter fasta file \(E\.g\. TruSeq3\-PE\.fa\) is required in the current working directory\.\n";
	print "\t\(Adapter fasta files can be downloaded from http\:\/\/github\.com\/usadellab\/Trimmomatic\)\n\n";
	
	print "Reference\:\n\tBolger\, A\. M\.\, Lohse\, M\.\, and Usadel\, B\. \(2014\)\. Trimmomatic\: A flexible trimmer for Illumina Sequence Data\. Bioinformatics\, btu170\.\n";

	exit;
}

##### II. Trimmomatic ################################################################################################################################

%hash = ();

opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;
@dir2 = @dir;

foreach $file2 (@dir2) { # Checking fastq files for R2 reads
	chomp $file2;
	if ($file2 =~ /_2\.fastq$/) {
		$file = $file2;
		$file =~ s/_2\.fastq$/_1\.fastq/;
		$hash{$file} = $file2;
	}
	if ($file2 eq "Fastq_output") {
		print "Directory for the output files already exists\n";
		exit;
	}
}

system ("mkdir \-m 777 Fastq_output");

foreach $file (@dir){
	chomp $file;
	if ($file =~ /_1\.fastq$/) {
		if (exists($hash{$file})) {
			$out1 = $file;
			$out1 =~ s/_1\.fastq$/_tm_1\.fastq/;
			$outunp1 = $out1;
			$outunp1 =~ s/_1\.fastq$/_unp_1\.fastq/;
			$out2 = $out1;
			$out2 =~ s/_tm_1\.fastq$/_tm_2\.fastq/;
			$outunp2 = $out2;
			$outunp2 =~ s/_2\.fastq$/_unp_2\.fastq/;
			if (length($adapterfa) > 0) {
				$sysin = "trimmomatic PE \-phred33 \-threads $threadnum $file $hash{$file} $out1 $outunp1 $out2 $outunp2 ILLUMINACLIP\:$adapterfa\:2\:30\:10 LEADING\:20 TRAILING\:20 SLIDINGWINDOW\:4\:15 MINLEN\:$minlen";
			} else {
				$sysin = "trimmomatic PE \-phred33 \-threads $threadnum $file $hash{$file} $out1 $outunp1 $out2 $outunp2 LEADING\:20 TRAILING\:20 SLIDINGWINDOW\:4\:15 MINLEN\:$minlen";
			}
			#print "$sysin\n";
			system($sysin);
			$sysin = "mv $out1 $outdir\/$out1";
			system($sysin);
			$sysin = "mv $out2 $outdir\/$out2";
			system($sysin);
			print "Move filtered fastq to the Fastq_output directory\n";
			unlink ($outunp1);
			unlink ($outunp2);
			print "Remove unpaired fastq\n\n";
		}
	}
}

exit;
