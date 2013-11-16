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
 
 
eval {
                CRI::SOAP->debug();
                $exportService = CRI::SOAP->new('http', $host, $port, '/solexa-ws/SolexaExportBeanWS', 'http://solexa.webservice.melab.cruk.org/');
};
if($@) {
                #print $@;
}




my $dir = getcwd;
open OUT, ">".$dir."/Projects_ActualLocations.txt";
open IN,"<".$dir."/Temp/Projects.txt";
my $line;
print "\n\n".$dir."/Temp/Projects.txt";
foreach $line (<IN>){
chomp($line);
print $line;
}
