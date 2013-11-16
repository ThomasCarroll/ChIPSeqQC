use strict;

my $L = $ARGV[0];
my $bed_filename = $ARGV[2];
my $input_filename = $ARGV[3];
my $DirOfChoice = $ARGV[4];
my $min_region_length = $ARGV[5];
my $width = $ARGV[1];
my $D = $ARGV[6];
#my $L = 200;
#my $D = 5;
#my $bed_filename = "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/bamFiles/SLX-4497.739.s_2.bwa.homo_sapiens_Processed.bed";
#my $min_region_length = 10000;
open OUT2, ">".$DirOfChoice."ZetaRegions.txt";
#my $input_filename = "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/bamFiles/SLX-4501.739.s_6.bwa.homo_sapiens_Processed.bed";
my @index_set = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");
foreach my $i (@index_set){
  my $chr = "chr".$i;
  my $outfile = $DirOfChoice."rho.".$chr.".txt";
  open IN, "<$input_filename";
    my %data = ();
    
    while (<IN>) {
	chomp;	
	my ($chrom, $pos, $endpos, $t1, $t2, $str, @rest) = split /\s+/, $_;
	if($chrom eq $chr){
	if ($str eq "+") {
		my $start = $pos-$width/2;
 		$start = 1 if ($start < 1);
		$data{$start} += 1;
		$data{$pos+$width} += -1;
	}
	elsif ($str eq "-") {
		my $start = $endpos - $width/2;
		 $start = 1 if ($start < 1);
		$data{$start} += 1;
		$data{$pos+$width}  += -1;
	}
	else {
		print "PROBLEM\t$str\n";
	}
	}
    }
    close IN;
 # Create a hash with all the read starts for the treatment
    open IN, "<$bed_filename";
    my %starts_treatment = ();

    while (<IN>) {
        chomp;
       my ($chrom, $pos, $endpos, $t1, $t2, $str, @rest) = split /\s+/, $_;
 	if($chrom eq $chr){
        	$pos = $endpos if ($str eq "-");
        	$starts_treatment{$pos} += 1;
	}
    }
    close IN;

   
    open OUT, ">$outfile";	
    my @sorted_pos = sort {$a <=> $b} (keys %data);
    my $height = 0;
    my $theta;
    my $begin_nuc=1;	
    foreach my $pos (@sorted_pos) {
	$height += $data{$pos};
	if($pos>=$width/2){
		$theta = $height * $L / $width;
	}else{
		$theta = $height * $L/($width/2 + $pos ) 
	}
	print OUT "$chr\t$pos\t$theta\n";
    }
    
    close OUT;
    open IN2, "$outfile";
	my @merged_regions = ();	
	my ($begin_nuc, $curr_nuc, $zz)= (1,1,0);
	my ($old_begin, $old_end, $old_rho, $end_nuc, $rho_avg);
	until(eof(IN2)){
	my($sum, $curr_rho) = (0,0);
	my $xx = 0;
	while($curr_nuc <= $begin_nuc + $min_region_length && $xx < 1){
		if(eof(IN2)){
			$xx = 500;
		}else{
			my($chromo,$next_nuc, $next_rho)= split /\s+/, <IN2>;
			$sum+=($next_nuc - $curr_nuc)*$curr_rho;			
			$curr_nuc = $next_nuc;
			$curr_rho = $next_rho;
		}
	}
	$rho_avg= $sum/($curr_nuc - $begin_nuc);
	my $diff = $D;
	if($rho_avg*.1 > $diff){
		$diff = $rho_avg*.1;
	}
	while(abs($curr_rho - $rho_avg) <= $diff && $xx<1){
		if(eof(IN2)){
			$xx = 500;
		}else{
			my($chromo,$next_nuc, $next_rho)= split /\s+/, <IN2>;
			$sum+=($next_nuc - $curr_nuc)*$curr_rho;			
			$curr_nuc = $next_nuc;
			$curr_rho = $next_rho;
		}
	}
	$rho_avg= $sum/($curr_nuc - $begin_nuc);
	$end_nuc = $curr_nuc -1;
	if($zz <1){
		($old_begin, $old_end, $old_rho)=($begin_nuc, $end_nuc, $rho_avg);
		$zz = 500;
	}else{
		if(int($old_rho+ .5) == int($rho_avg + .5)){
			$old_end = $end_nuc;
		}else{
			my $value =int($old_rho + .5)+1;
			push @merged_regions, "$old_begin\t$old_end\t$value";
			($old_begin, $old_end, $old_rho)=($begin_nuc, $end_nuc, $rho_avg);
		}
	}
	$begin_nuc = $curr_nuc; 
	}
	my $value =int($old_rho+ .5)+1;
       	push @merged_regions, "$old_begin\t$end_nuc\t$value";

    my @treatment_positions = sort {$a <=> $b} keys %starts_treatment;   # Sorted list of treatment start positions
    my $no_treatment_positions = @treatment_positions;
foreach (@merged_regions){
	my ($start, $stop,$theta) = split /\s+/, $_; 	
	
        my $tags_treatment          = 0;

        # Count the number of treatment tags in the region

        my $index = int ($no_treatment_positions / 2);
        my $shift = int ($no_treatment_positions / 4);
        my $flag = 0;

        while ($flag == 0) {
            my $pos      = $treatment_positions[$index];
            my $prec_pos = $treatment_positions[$index - 1];

            if ($prec_pos < $start - $L && $pos >= $start - $L) {
                $flag = 1;
            }
            elsif ($index == 1 || $index == $no_treatment_positions - 1) {
                $flag = 1;
            }
            elsif ($pos < $start - $L) {
                $index += $shift;
                $shift = max(1, int ($shift / 2));
            }
            elsif ($prec_pos >= $start - $L) {
                $index -= $shift;
                $shift = max(1, int ($shift / 2));
            }
        }

        while ($treatment_positions[$index] <= $stop && $index < $no_treatment_positions) {
            $tags_treatment += $starts_treatment{$treatment_positions[$index]};
            $index++;
        }
 print OUT2 "$chr\t$start\t$stop\t$tags_treatment\t$theta\n";
}
	system("rm $outfile"); 
    print "$chr\n";
    
}
close OUT2;
sub max {
    my ($a, $b) = @_;
    if ($a > $b) {
        return $a;
    }
 else {
        return $b;
    }
}
