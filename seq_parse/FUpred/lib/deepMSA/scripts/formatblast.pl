#!/usr/bin/perl
use Math::BigFloat;
use Math::Trig;

#$blastdir="/nfs/amino-home/zhengwei/Blast/blast-2.2.26/bin/";
#$db="/nfs/amino-library/nr/nr";
############# example: ./formatblast.pl d1sdwa2.fasta blast.out msa.txt 100 0.001
$fastafile=$ARGV[0];
$blastfile=$ARGV[1];
$outfile=$ARGV[2];
$idcut=$ARGV[3];
$idcut=~/(\d+.*\d)/;
$idcut=$1;
print "$idcut\n";
$evaluecut=$ARGV[4];
$evaluecut=~/(\d+.*\d)/;
$evaluecut=$1;
print "$evaluecut\n";
################ make fasta sequence file #################
@seqtxts=`cat $fastafile`;
$sequence="";
foreach $seqtxt(@seqtxts){
    goto pos6 if($seqtxt=~/\>/);
    $seqtxt=~s/\s//mg;
    $seqtxt=~s/\n//mg;
    $sequence=$sequence.$seqtxt;
  pos6:;
}
$Lch=length $sequence;

########### run psi-blast #######################
#printf "running Psi-blast .....\n";
#`$blastdir/blastpgp  -b 1000 -j 3 -h 0.001 -d $db -i protein.seq -C psitmp.chk > blast.out`;

#### record multiple sequence alignment $am{i_seq,i_pos} -------->
open(blast,$blastfile);
while($line=<blast>){
    if($line=~/Results from round\s+(\d+)/){
        $ROUND=$1;
    }
}
seek(blast,0,0);
$it=0;
print "$ROUND\n";
while($line=<blast>){
    if($line=~/Results from round\s+$ROUND/ || $ROUND eq ''){
        while($line=<blast>){
	    if($line=~/>(\S+)/){
		$name{$it}=$1;
	    }
            if($line=~/Expect =\s*(\S+)\,/){
		if($name{$it}  eq ''){
		    $name{$it}=$name{$it-1};
		}
                $ev=$1;
		if($ev=~/^e\S+/){
			$ev1="1".$ev;
		}
		else{
			$ev1=$ev;
		}
		$evalue{$it}=$ev1;
		$ev1 = new Math::BigFloat $ev1;

		#print "$name{$it}\n";
		#print "$ev1\n";

                <blast>=~/Identities =\s*\S+\s+\((\S+)\%/;
                <blast>;
                $id=$1;
                if($ev1<$evaluecut && $id < $idcut){
		    $it++; #number of aligned sequences
                    while($line=<blast>){
                        if($line=~/Query\:\s*(\d+)\s+(\S+)\s+(\S+)/){
                            $i1=$1;
                            $seq1=$2;
                            <blast>;
                            <blast>=~/Sbjct\:\s*(\S+)\s+(\S+)\s+(\S+)/;
                            $seq2=$2;
                            <blast>;
                            ###
                            $L=length $seq1;
                            $m1=$i1-1;
                            for($i=1;$i<=$L;$i++){
                                $q1=substr($seq1,$i-1,1);
                                $q2=substr($seq2,$i-1,1);
                                $m1++ if($q1 ne '-');
                                if($q1 ne '-'){
				    $am{$it,$m1}=$q2;
                                }
                            }
                            ###
                        }else{
                            goto pos1;
                        }
                    }
                }
              pos1:;
            }
        }
    }
}
close(blast);

open(msa,">$outfile");
printf msa ">$fastafile\n";
printf msa "$sequence\n";
for($j=1;$j<=$it;$j++){
    printf msa ">$name{$j-1}#psi-ev:$evalue{$j-1}\n";
    for($i=1;$i<=$Lch;$i++){
	if($am{$j,$i} ne ''){
	   printf msa "$am{$j,$i}";
	}
        else{
	   printf msa "-";
	}
    }
    printf msa "\n";
}
close(msa);



