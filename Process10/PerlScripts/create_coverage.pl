use strict;

my $L = $ARGV[0];
my $bed_filename = $ARGV[1];
my $DirOfChoice = $ARGV[3];
my $Genome = $ARGV[2];
my @index_set;
my $Indexsetfile = $ARGV[4];

open (IN,$Indexsetfile);
my $count = 0;
my %data = ();
#print $count;
while (<IN>) {
	chomp;	
	my ($chrom, $len) = split /\s+/, $_;
	#print $chrom;
	$index_set[$count] = $chrom;
	$count++;
	#print $count;
	}
close(IN); 
#print @index_set;
foreach my $i (@index_set){
  my $chr = $i;
  my $outfile = $DirOfChoice."Coverage.".$chr.".sgr";
  open OUT, ">$outfile";
  open IN, "$bed_filename";
  my %data = ();
    
    while (<IN>) {
	chomp;
	my ($chrom, $pos, $endpos, $t1, $t2, $str, @rest) = split /\s+/, $_; 
	if($chrom eq $chr){
	if ($str eq "+") {
		$data{$pos} += 1;
		$data{$pos+$L} += -1;
	}
	elsif ($str eq "-") {
		my $start = $endpos - $L;
		$start = 1 if ($start < 1);
		my $stop  = $endpos;
		$data{$start} += 1;
		$data{$stop}  += -1;
	}
	else {
		print "PROBLEM\n";
	}
	}
    }
    
    close IN;
   


	
    my @sorted_pos = sort {$a <=> $b} (keys %data);
    my $height = 0;
	
    foreach my $pos (@sorted_pos) {
	$height += $data{$pos};
	print OUT "$chr\t$pos\t$height\n";
    }
    
    close OUT;
    print "$chr\n";
    
}
