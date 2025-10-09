
#!/usr/bin/perl

# generates and runs three sets of R scripts for input to spatgen

$nl = 1000;
$np = 10;
$ng = 10;
$nt = $np*$ng;
$nrep = 10;
$ne = 1000;
$N = 100;

### constant selection sims ###
for($r=0; $r<$nrep; $r++){
    open(R, "> source.R") or die "failed to write source.R\n";
    print R "linearCurve<-function(A=NA,B=NA,x=seq(-1,1,0.01)){\n";
    print R "y<-A + B*x\n}\n";
    print R "evar<-rnorm(100)\n";
    print R "p<-runif($nl,0.05,0.95)\n";
    print R "P<-vector(\"list\",$nl)\n";
    print R "for(i in 1:$nl){P[[i]]<-matrix(p[i],nrow=$np,ncol=$ng)}\n";
    print R "A<-P\n";
    print R "for(k in 1:($ng-1)){\nfor(i in 1:$nl){\nfor(j in 1:$np){\n";
    print R "p<-P[[i]][j,k]\n";
    print R "A[[i]][j,k]<-rbinom(1,size=2*$N,p)\n";
    print R "if(i==1){\n";## selected locus, fitness based on environment
    print R "X<-p^2\nY2<-2*p*(1-p)\nZ<-(1-p)^2\n";
    print R "x<-k + (j-1) * $ng\n";
    print R "w<-rep(NA,3)\n";
    print R "w[1]<-linearCurve(A=0.9,B=0.1,evar[x])\n";
    print R "w[2]<-linearCurve(A=0.9,B=0,evar[x])\n";
    print R "w[3]<-linearCurve(A=0.9,B=-0.1,evar[x])\n";
    print R "wbar<-X*w[1] + Y2 * w[2] + Z * w[3]\n";
    print R "p<-p + (X*w[1] + 0.5 * Y2 * w[2] - p * wbar)/wbar\n}\n";
    print R "p<-rbinom(1,size=2*$ne,p)/(2*$ne)\n";
    print R "p[p<0.001]<-0.001\np[p>0.999]<-0.999\n";
    print R "P[[i]][j,k+1]<-p\n";
    print R "}\n}\n}\n";
    print R "for(i in 1:$nl){\nfor(j in 1:$np){\n";
    print R "k=$ng\n";
    print R "p<-P[[i]][j,k]\n";
    print R "A[[i]][j,k]<-rbinom(1,size=2*$N,p)\n";
    print R "}\n}\n";
    print R "G<-matrix(NA,nrow=$nl,ncol=$nt)\n";
    print R "X<-matrix(NA,nrow=$nl,ncol=$nt)\n";
    print R "for(i in 1:$nl){G[i,]<-as.vector(t(A[[i]]))}\n";
    print R "for(i in 1:$nl){X[i,]<-as.vector(t(P[[i]]))}\n";
    print R "write.table(round(evar,5),\"rsimFluctLE"."$r".".txt\",row.names=F,quote=F,col.names=F)\n";
    print R "write.table(G,\"rsimFluctLG"."$r".".txt\",row.names=F,col.names=F,quote=F)\n";
    print R "write.table(round(X,4),\"rsimFluctLP"."$r".".txt\",row.names=F,col.names=F,quote=F)\n";
    close(R);
    system "R CMD BATCH source.R\n";

## convert allele frequency file from R to spatgen

    $in = "rsimFluctLG"."$r".".txt";
    open (IN, $in) or die "failed to read\n";
    open (OUT, "> mod_$in") or die "failed to write\n";
    print OUT "$np $ng $nl\n";;
    while(<IN>){
	chomp;
	@line = split(" ",$_);
	@a = ();
	foreach $p (@line){
	    $x = 2*$N - $p;
	    push (@a, "$p,$x");
	}
	$a = join(" ",@a);
	print OUT "$a\n";
    }
    close(IN);
    close(OUT);
}

