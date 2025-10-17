#!/usr/bin/perl
use File::Path qw(make_path);

my $qbin = $ENV{QNEMO_BIN} // 'quantinemo';   # or 'quantiNemo2' / 'quantiNemo'

# wrapper script for quantiNemo
# generates selection file, runs quantiNemo, and extracts necessary output

##### parameters ########
$npops = 10;
$ngens = 20;
$nloci = 1000;
$nqloci = 10;
$a = 0.2; # additive effect
$spat_mn_sel = 0;
$spat_sd_sel = 1;
$temp_sd_sel = 1;
$w = 1.5; # selection intensity
## 0-4 have 2
#########################


$nreps = 10;

for($i=9; $i<$nreps; $i++){
    # ## use R to generate matrix optimal phenotypes
    open(R, "> source$i.R") or die "failed to write the R file\n";
    print R "mns<-rnorm(n=$npops,mean=$spat_mn_sel,sd=$spat_sd_sel)\n";
    print R "env<-matrix(NA,nrow=$ngens,ncol=$npops)\n";
    print R "for(j in 1:$npops){env[,j]<-rnorm(n=$ngens,mean=mns[j],sd=$temp_sd_sel)}\n";
    print R "write.table(round(env,5),\"file$i.txt\",row.names=FALSE,col.names=FALSE,quote=FALSE)\n";
    close(R);
    system "R CMD BATCH source$i.R\n";

    ## convert this into the file format for quantiNemo
    open (SFILE, "file$i.txt") or die "failed to read output from R\n";
    open (SOUT, "> selection_file$i.txt") or die "failed to write selection_file.txt\n";
    $t = 1;
    print SOUT "(";
    while (<SFILE>){
	chomp;
	if ($t < $ngens){
	    print SOUT "$t {$_},\n";
	}
	else{
	    print SOUT "$t {$_})\n";
	}
	$t++;
    }


    ## write the quanti_allelic_file
    open (QINI, "> q_ini$i.txt") or die "failed to write QINI\n";
    print QINI "[FILE_INFO]{\n col_locus 1\n col_allele 2\n col_allelic_value 3\n col_ini_freq 4\n}\n";
    for($l=1; $l<=$nqloci; $l++){
	$p = 0.1 + rand(0.8);
	$q = 1 - $p;
	print QINI "$l 1 -$a $p\n$l 2 $a $q\n";
    }
    close (QINI);
    
    ## write the ntrl_allelic_file
    open (NINI, "> ntrl_ini$i.txt") or die "failed to write NINI\n";
    print NINI "[FILE_INFO]{\n col_locus 1\n col_allele 2\n col_ini_freq 3\n}\n";
    for($l=1; $l<=$nloci; $l++){
        $p = 0.05 + rand(0.9);
        $q = 1 - $p;
        #print NINI "$l 1 $p\n$l 2 $q\n";

        print NINI "$l 1 $p";
        for($j=2; $j<99999; $j++){ print NINI " 0";}
        print NINI "\n";

        print NINI "$l 2 $q";
        for($j=2; $j<99999; $j++){ print NINI " 0";}
        print NINI "\n";
    }
    close (NINI);

    ## write the quantiNemo file and run it
    open (SIM, "> fluct$i.ini") or die "could not write fluct.ini\n";

    print SIM "logfile_type              2\n";
    print SIM "folder		     simout$i\n";
    print SIM "overwrite		     1\n"; ## overwrite automatically
    print SIM "generations               $ngens\n";
    print SIM "genetic_map_output	     1\n";

    # metapopulation
    print SIM "patch_number              $npops\n";
    print SIM "patch_capacity            1000\n";
    print SIM "dispersal_rate            0.001\n";

    # mating
    print SIM "breed_model               0\n";
    print SIM "mating_system             0\n";

    # selection
    print SIM "patch_stab_sel_optima     \$selection_file$i.txt\n";
    print SIM "patch_stab_sel_intensity  $w\n";

    # quantitative trait
    print SIM "quanti_loci               $nqloci\n";
    print SIM "quanti_loci_positions     {";
    foreach $chr (1..5) {
	print SIM "{$chr".": 25.5 75.5}\n";
    }
    print SIM "}\n";
    print SIM "quanti_all                2\n";
    print SIM "quanti_mutation_model     0\n";
    print SIM "quanti_mutation_rate      1e-4\n";
    print SIM "quanti_allelic_file       q_ini$i.txt\n";
    print SIM "quanti_save_genotype      1\n";
    print SIM "quanti_genot_logtime      (1 20, 10 1)\n";
    print SIM "quanti_genot_dir          quanti_genotype\n";

    print SIM "quanti_save_phenotype     1\n";
    print SIM "quanti_phenot_logtime     (1 20, 10 1)\n";
    print SIM "quanti_phenot_dir         quanti_phenotype\n";
    print SIM "quanti_output             1\n";

    print SIM "quanti_save_geno_value    1\n";
    print SIM "quanti_geno_value_logtime (1 20, 10 1)\n";
    print SIM "quanti_geno_value_dir     quanti_genoval\n";


    # neutral marker
    print SIM "ntrl_loci                 $nloci\n";
    print SIM "ntrl_loci_positions     {";
    foreach $chr (1..10) {
	print SIM "{$chr".":";
	foreach $pos (1..100) {
	    print SIM " $pos";
	}
	print SIM "}\n";
    }
    print SIM "}\n";

    print SIM "ntrl_all                  2\n";
    print SIM "ntrl_mutation_model       0\n";
    print SIM "ntrl_mutation_rate        1e-4\n";
    print SIM "ntrl_allelic_file         ntrl_ini$i.txt\n";

    print SIM "ntrl_save_genotype        1\n";
    print SIM "ntrl_genot_logtime         (1 20, 10 1)\n";
    print SIM "ntrl_genot_dir            ntrl_genotype\n";

    # statistics
    print SIM "stat                      {VwW meanW_p}\n";
    print SIM "stat_save                 1\n";
    print SIM "stat_log_time             (1 20, 10 1)\n";
    print SIM "stat_dir                  stats\n";

    close(SIM);
    #system "quantiNemo fluct$i.ini\n";
    my $rc = system("$qbin fluct$i.ini");
    die "quantiNemo run failed for replicate $i (exit ".($? >> 8).")\n" if $rc != 0;

    ## grab genotypic values
    open (OUT, "> simout$i"."/breedingval.txt") or die "failed to write OUT\n";
    foreach $gen (10..20){
	##  file
	if($gen == 20){
	    $file = "simout$i/quanti_genoval/simulation_g20.gen";
	}
	else{
	    $file = "simout$i/quanti_genoval/simulation_g$gen".".gen";
	}
	open(IN, $file) or die "could not open the genotype file: $file\n";
	<IN>;
	<IN>;
	while(<IN>){ 
	    m/^(\d+)\s+(\S+)/;
	    $pop = $1;
	    $val = $2;
	    push (@{$gval{$pop}{$gen}}, $val); 
	    print OUT "$pop $gen $val\n";
	}
	close(IN);
    }
    close(OUT);

    ## grab optimal phenotype for each generation
    open(IN, "selection_file$i.txt") or die "couldn't read the selection file\n";
    while (<IN>){
	chomp;
	if (m/^(1\d)/ | m/^(20)/){
	    $gen = $1;
	    s/[{,}]//g;
	    @data = split(" ",$_);
	    shift(@data);
	    $pop = 1;
	    foreach $opt (@data){
		$selopt{$pop}{$gen} = $opt;
		$pop += 1;
	    }
	}
    }
    close (IN);

    ## calculate each individual's expected fitness
    open (OUT, "> simout$i"."/fitness.txt") or die "failed to write OUT\n";
    foreach $pop (sort {$a<=>$b} keys %gval){
	foreach $gen (sort {$a<=>$b} keys %{$gval{$pop}}){
	    foreach $ind (@{$gval{$pop}{$gen}}){
		$val = -1 * ((($ind - $selopt{$pop}{$gen}) ** 2)/(2 * $w ** 2));
		$val = exp($val);
	        push (@{$fitness{$pop}{$gen}},$val);
		print OUT "$pop $gen $val\n";
	    }
	}
    }
    close (OUT);

    ## calculate sample allele frequencies for each population,
    ## generation and locus; also calculate selection experienced by
    ## each locus, i.e. difference in marginal fitness of alternative
    ## alleles
    open (OUT, "> simout$i"."/quantiGenotypes.txt") or die "failed to write OUT\n";
    foreach $gen (10..20){
	## quanti_genotype file
	if($gen == 20){
	    $file = "simout$i/quanti_genotype/simulation_g20.dat";
	}
	else{
	    $file = "simout$i/quanti_genotype/simulation_g$gen".".dat";
	}
	open(IN, $file) or die "could not open the genotype file: $file\n";
	<IN>;
	while(<IN>){
	    chomp;
	    if (m/^\d+/){## genotype lines
		@data = split(" ",$_);
		$pop = shift(@data);
		$locus = 0;
		print OUT ("$pop $gen");
		foreach $G (@data){
		    @G = split("",$G);
		    $g = $G[0] + $G[1] - 2;
		    $allfreq{$locus}{$pop}{$gen} += $g;
		    print OUT " $g";
		    $cnt{$locus}{$pop}{$gen} += 2;
		    $locus++;
		}
		print OUT "\n";
	    }
	}
	close(IN);
    }
    close (OUT);

    open (OUT, "> simout$i"."/QuantiAlleleFreqsSim.txt") or die "failed to write OUT\n";
    foreach $locus (sort {$a<=>$b} keys %cnt){
	foreach $pop (sort {$a<=>$b} keys %{$cnt{$locus}}){
	    foreach $gen (sort {$a<=>$b} keys %{$cnt{$locus}{$pop}}){
		#print "l: $locus; p: $pop; g: $gen\n";
		print OUT "$allfreq{$locus}{$pop}{$gen} ";
	    }
	}
	print OUT "\n";
    }
    close (OUT);
    
    ## allele frequency change by environment
    
    open (OUT, "> simout$i"."/quantiDp.txt") or die "could not write the outfile\n";
    $locus = 1;
    foreach $pop (sort {$a<=>$b} keys %{$cnt{$locus}}){
	foreach $gen (sort {$a<=>$b} keys %{$cnt{$locus}{$pop}}){
	    unless ($gen == 20){
	        $p = $pop;
		$p =~ s/^0//;
		$env = $selopt{$p}{$gen};
		#print "$gen $env\n";
		$g1 = $gen;
		#$g1 =~ s/^0//;
		print OUT "$pop $gen $env";
		foreach $locus (sort {$a<=>$b} keys %cnt){
		    $g2 = $g1 + 1;

		    $dp = $allfreq{$locus}{$pop}{$g2}/2000 - $allfreq{$locus}{$pop}{$g1}/2000;
		    print OUT " $dp";
		}
		print OUT "\n";
	    }
	}
    }

    ## selection by each locus
    open(R, "> simout$i"."/source.R") or die "failed to write an R file\n";
    print R "w<-as.matrix(read.table(\"simout$i"."/fitness.txt\",header=FALSE))\n";
    print R "G<-as.matrix(read.table(\"simout$i"."/quantiGenotypes.txt\",header=FALSE))\n";
    print R "s<-matrix(NA,nrow=$npops * 11,ncol=2 + $nqloci)\n";
    print R "n<-1\n";
    print R "for(j in 1:$npops){\nfor(t in 10:20){\n";
    print R "x<-(w[,1] == j & w[,2] == t)\n";
    print R "y<-(G[,1] == j & G[,2] == t)\n";
    print R "for(l in 1:$nqloci){\n";
    print R "if(length(unique(G[y,l+2])) > 1){\n";
    print R "out<-lm(w[x,3] ~ G[y,l+2])\n";
    print R "s[n,1]<-j\ns[n,2]<-t\ns[n,l+2]<-out\$coefficients[2]\n";
    print R "}\n";
    print R "else{s[n,1]<-j\ns[n,2]<-t\ns[n,l+2]<-NA\n}\n}\n";
    print R "n<-n+1\n";
    print R "}\n}\n";
    print R "write.table(round(s,4),\"simout$i"."/quantiLocusSelection.txt\",row.names=FALSE,col.names=F,quote=FALSE)\n";
    
    close(R);

    system "R CMD BATCH simout$i"."/source.R\n";

    %allfreq = ();
    %cnt = ();


    open (OUT, "> simout$i"."/ntrlGenotypes.txt") or die "failed to write OUT\n";
    ## repeat for neutral loci 
    foreach $gen (10..20){
	## neutral_genotype file
	if($gen == 20){
	    $file = "simout$i/ntrl_genotype/simulation_g20.dat";
	}
	else{
	    $file = "simout$i/ntrl_genotype/simulation_g$gen".".dat";
	}
	open(IN, $file) or die "could not open the genotype file: $file\n";
	<IN>;
	while(<IN>){
	    chomp;
	    if (m/^\d+/){## genotype lines
		@data = split(" ",$_);
		$pop = shift(@data);
		$locus = 0;
		print OUT ("$pop $gen");
		foreach $G (@data){
		    @G = split("",$G);
		    $g = $G[0] + $G[1] - 2;
		    $allfreq{$locus}{$pop}{$gen} += $g;
		    print OUT " $g";
		    $cnt{$locus}{$pop}{$gen} += 2;
		    $locus++;
		}
		print OUT "\n";
	    }
	}
	close(IN);
    }
    close (OUT);

    

    open (OUT, "> simout$i"."/NtrlAlleleFreqsSim.txt") or die "failed to write OUT\n";
    foreach $locus (sort {$a<=>$b} keys %cnt){
	foreach $pop (sort {$a<=>$b} keys %{$cnt{$locus}}){
	    foreach $gen (sort {$a<=>$b} keys %{$cnt{$locus}{$pop}}){
		#print "l: $locus; p: $pop; g: $gen\n";
		print OUT "$allfreq{$locus}{$pop}{$gen} ";
	    }
	}
	print OUT "\n";
    }
    close (OUT);


    ## allele frequency change by environment
    
    open (OUT, "> simout$i"."/ntrlDp.txt") or die "could not write the outfile\n";
    $locus = 1;
    foreach $pop (sort {$a<=>$b} keys %{$cnt{$locus}}){
	foreach $gen (sort {$a<=>$b} keys %{$cnt{$locus}{$pop}}){
	    unless ($gen == 20){
	        $p = $pop;
		$p =~ s/^0//;
		$env = $selopt{$p}{$gen};
		$g1 = $gen;
		#$g1 =~ s/^0//;
		print OUT "$pop $gen $env";
		foreach $locus (sort {$a<=>$b} keys %cnt){
		    $g2 = $g1 + 1;

		    $dp = $allfreq{$locus}{$pop}{$g2}/2000 - $allfreq{$locus}{$pop}{$g1}/2000;
		    print OUT " $dp";
		}
		print OUT "\n";
	    }
	}
    }


    ## selection by each locus
    open(R, "> simout$i"."/source.R") or die "failed to write an R file\n";
    print R "w<-as.matrix(read.table(\"simout$i"."/fitness.txt\",header=FALSE))\n";
    print R "G<-as.matrix(read.table(\"simout$i"."/ntrlGenotypes.txt\",header=FALSE))\n";
    print R "s<-matrix(NA,nrow=$npops * 11,ncol=2 + $nloci)\n";
    print R "n<-1\n";
    print R "for(j in 1:$npops){\nfor(t in 10:20){\n";
    print R "x<-(w[,1] == j & w[,2] == t)\n";
    print R "y<-(G[,1] == j & G[,2] == t)\n";
    print R "for(l in 1:$nloci){\n";
    print R "if(length(unique(G[y,l+2]))>1){\nout<-lm(w[x,3] ~ G[y,l+2])\n";
    print R "s[n,1]<-j\ns[n,2]<-t\ns[n,l+2]<-out\$coefficients[2]\n";
    print R "}\n";
    print R "else{s[n,1]<-j\ns[n,2]<-t\ns[n,l+2]<-NA\n}\n}\n";
    print R "n<-n+1\n";
    print R "}\n}\n";
    print R "write.table(round(s,4),\"simout$i"."/ntrlLocusSelection.txt\",row.names=FALSE,col.names=F,quote=FALSE)\n";
    
    close(R);

    system "R CMD BATCH simout$i"."/source.R\n";

}

