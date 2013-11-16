#!/usr/bin/perl -w
 
use strict;
use CRI::FileLocation;
use CRI::SOAP;
use constant {TRUE => 1, FALSE => 0};
use File::Basename; 
use Cwd;
use File::Spec;

my ($exportService);
my $host = 'uk-cri-ldmz02.crnet.org';
my $port = 80;
 


my $dir = getcwd;

open IN,"<".$dir."/SLXToAlign.txt";
my $line;
foreach $line (<IN>){
	print $line;
	chomp($line);
	print $line;
}
