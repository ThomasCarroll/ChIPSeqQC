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
#print $dir;
my $FileOut=join($dir."ActualLocations.txt");
open OUT, ">".$dir."/ActualLocations.txt";
#print $FileOut;
my $Project = basename($dir);
#my $lanes = $exportService->call('getAllSampleLanesForId','SLX-2167');
print $Project;
my $lanes = $exportService->call('getLanesFromAnnotation','Project',$Project);
foreach my $lane (@$lanes) {
		#print $lane;
                my $locations = $exportService->call('getFileLocations', $lane->{'sampleProcess_id'}, 'FILE', 'FASTQ');
		foreach my $location (@$locations) {
                                printHash($location);
                }

}
close OUT; 
 
sub printHash {
                my ($hash) = @_;
                while ( (my $k,my $v) = each %$hash) {
                                print OUT " $k => $v\n";
					print  " $k => $v\n";
                }
                print OUT "\n";
}

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
