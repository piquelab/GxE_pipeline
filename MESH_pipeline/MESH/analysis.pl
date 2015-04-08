# input mesh.in from Chris
$ARGV[0] =~/\/(\S+)\.txt/;

$prefix = $1;

## modify the format
## not needed if chris fixes the R script
open OUT, ">$prefix.mesh.in";
while(<>){
    next if $_ !~ /\d/;
    
    my @data = split/\s+/, $_;
    shift @data until @data[0]=~/^\S/;
    
    next if $data[0] eq '.';

    $data[3] = 1 if($data[3] eq "NaN");
    print OUT "$data[0]\_$data[0]  $data[1] $data[2] $data[3]\n";
}

## #######################################

#### change $prefix.mesh.in to $ARGV[0], once the above is fixed

`mesh -d $prefix.mesh.in -g grid -min_info -prep_hm > $prefix.hm_config.in`;
`mesh -d $prefix.mesh.in -g grid.het -min_info -prep_hm > $prefix.temp`;

#`mesh -d $ARGV[0] -g grid -min_info -prep_hm > $prefix.hm_config.in`;
#`mesh -d $ARGV[0] -g grid.het -min_info -prep_hm > $prefix.temp`;



open FILE, "$prefix.temp";
open OUT, ">$prefix.hm_het.in";
while(<FILE>){
    
    next if $_ !~ /^\s*\S+\s+(\d)\s+/;
    next if $1 != 3;
    print OUT $_;
}


unlink "$prefix.temp";

print STDERR "running config MCMC ...\n";
`hm_mh -d $prefix.hm_config.in -g 3 -s 3 -b 10000 -r 10000 > $prefix.hm_config.out 2>$prefix.hm_config.est`;

print STDERR "running het MCMC ...\n";
`hm_mh -d $prefix.hm_het.in -g 33 -s 1 -b 10000 -r 10000 > $prefix.hm_het.out 2>$prefix.het.est`;

open FILE, "$prefix.het.est";
open OUT, ">$prefix.hm_het.est";
while(<FILE>){
    next if $_ !~ /^\s*grid\:/;
    s/(\[\S+\s+\S+\])/ /g;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    shift @data;
    print OUT "posterior correlation: ";
    for ($i=0;$i<=$#data;$i+=6){
	printf OUT "%7.3e  ",$data[$i]+$data[$i+2]+$data[$i+4];
    }
    print OUT "\n";
}

#unlink "$preix.het.est", "$prefix.mesh.in", "$prefix.hm_het.in", "$prefix.hm_config.in";
