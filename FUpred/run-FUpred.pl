#!/usr/bin/perl
use Math::Trig;
use File::Basename;
use Cwd 'abs_path';

######## set your variables here
#$msa_hhblitsdb="/nfs/amino-home/zhng/local_library/uniclust30_2017_04/uniclust30_2017_04";
$msa_hhblitsdb="/public/dbu/data/af2_db/uniclust30/uniclust30_2018_08/uniclust30_2018_08";
#$msa_jackhmmerdb="/nfs/amino-home/zhng/local_library/uniref90_04042018/uniref90.fasta";
$msa_jackhmmerdb="/public/dbu/data/fupred_db/uniref90_db/uniref90.fasta";

$msa_hmmsearchdb="/nfs/amino-home/zhng/local_library/metaclust/metaclust.fasta";
#$python="/nfs/amino-home/zhng/local_library/anaconda3/bin/python";
$python="/public/dbu/local/anaconda3/envs/pytorch/bin/python";
#$python2="/nfs/amino-home/zhng/local_library/anaconda2/bin/python";
$python2="/public/dbu/local/anaconda3/envs/py27/bin/python";
$cpunum=4;

#################################
##### report node -------->
`hostname`=~/(\S+)/;
$node=$1;
printf "hostname: $node\n";
$time=`date`;
printf "starting time: $time";
$pwd=`pwd`;
printf "pwd: $pwd";
#^^^^^^^^^^^^^^^^^^^^^^^^^^

%ts=(
     'GLY'=>'G',
     'ALA'=>'A',
     'VAL'=>'V',
     'LEU'=>'L',
     'ILE'=>'I',
     'SER'=>'S',
     'THR'=>'T',
     'CYS'=>'C',
     'MET'=>'M',
     'PRO'=>'P',
     'ASP'=>'D',
     'ASN'=>'N',
     'GLU'=>'E',
     'GLN'=>'Q',
     'LYS'=>'K',
     'ARG'=>'R',
     'HIS'=>'H',
     'PHE'=>'F',
     'TYR'=>'Y',
     'TRP'=>'W',

     'ASX'=>'B',
     'GLX'=>'Z',
     'UNK'=>'X',

     'G'=>'GLY',
     'A'=>'ALA',
     'V'=>'VAL',
     'L'=>'LEU',
     'I'=>'ILE',
     'S'=>'SER',
     'T'=>'THR',
     'C'=>'CYS',
     'M'=>'MET',
     'P'=>'PRO',
     'D'=>'ASP',
     'N'=>'ASN',
     'E'=>'GLU',
     'Q'=>'GLN',
     'K'=>'LYS',
     'R'=>'ARG',
     'H'=>'HIS',
     'F'=>'PHE',
     'Y'=>'TYR',
     'W'=>'TRP',

     'a'=>'CYS',
     'b'=>'CYS',
     'c'=>'CYS',
     'd'=>'CYS',
     'e'=>'CYS',
     'f'=>'CYS',
     'g'=>'CYS',
     'h'=>'CYS',
     'i'=>'CYS',
     'j'=>'CYS',
     'k'=>'CYS',
     'l'=>'CYS',
     'm'=>'CYS',
     'n'=>'CYS',
     'o'=>'CYS',
     'p'=>'CYS',
     'q'=>'CYS',
     'r'=>'CYS',
     's'=>'CYS',
     't'=>'CYS',
     'u'=>'CYS',
     'v'=>'CYS',
     'w'=>'CYS',
     'x'=>'CYS',
     'y'=>'CYS',
     'z'=>'CYS',

     'B'=>'ASX',
     'Z'=>'GLX',
     'X'=>'CYS',
    );

@AA=qw(
       C
       M
       F
       I
       L
       V
       W
       Y
       A
       G
       T
       S
       Q
       N
       E
       D
       H
       R
       K
       P
       );

################# directories #############################
$infile=abs_path($ARGV[0]);
#print "$infile\n";
#$datadir=$ENV{'PWD'};
#$datadir=abs_path($ARGV[0]);
#$infile="$datadir/seq.fasta";
$scriptsdir=dirname(abs_path(__FILE__));

$bindir=abs_path("$scriptsdir/./bin");
$fudir=abs_path("$bindir/../");
#print "$bindir\n";
$method="respre";
$work_dir="/tmp/$ENV{'USER'}/FUD_".&getRandomString(10);

$psipreddir="$fudir/lib/psipred-4.01";
$ENV{'HHLIB'}="$fudir/lib/deepMSA";
$hhsuite=$ENV{'HHLIB'};
#print "$hhsuite\n";
################ working directory ########################
`/bin/mkdir -p $work_dir`;
chdir "$work_dir";
`/bin/rm -rf $work_dir/*`;
################ make fasta sequence file #################
@seqtxts=`cat $infile`;
$sequence="";
foreach $seqtxt(@seqtxts){
    goto pos1 if($seqtxt=~/\>/);
    $seqtxt=~s/\s//mg;
    $seqtxt=~s/\n//mg;
    $sequence=$sequence.$seqtxt;
  pos1:;
}
$Lch=length $sequence;
open(seq,">protein.seq");
printf seq ">protein\n";
for($i=1;$i<=$Lch;$i++)
{
    $a=substr($sequence,$i-1,1);
    foreach $A(@AA)
    {
	goto pos2 if($a eq $A);
    }
    $a="G";
  pos2:;
    printf seq "$a";
    $seqQ{$i}=$a;   #only for check
    $log{$i,$seqQ{$i}}++;
    if($i==int($i/60)*60)
    {
	printf seq "\n";
    }
}
printf seq "\n";
close(seq);

########### get mtx and prf from deepmsa ######################

system("$python2 $hhsuite/scripts/build_MSA.py -ncpu=$cpunum protein.seq -hhblitsdb=$msa_hhblitsdb -jackhmmerdb=$msa_jackhmmerdb  -hmmsearchdb=$msa_hmmsearchdb -tmpdir=$work_dir/deepMSA");

`cp protein.aln deepmsa_protein.aln`;
`cp protein.a3m deepmsa_protein.a3m`;
#`cp protein.aln $datadir/deepmsa.aln`;
#`cp protein.a3m $datadir/deepmsa.a3m`;

`$ENV{'HHLIB'}/scripts/a3m2psiblast_v2.pl deepmsa_protein.a3m psitmp -a3m -neff 7`;   ## for ss

`cp psitmp.mtx  protein.mtx`;
`cp psitmp.pssm protein.pssm`;

open(msa,"deepmsa_protein.aln");
$it=0;
my $queryseq=<msa>;
$L=length $queryseq;
$it++;
while($line=<msa>)
{
    for($i=1;$i<=$L;$i++)
    {
        $am{$it,$i}=substr($line,$i-1,1);
    }
    $it++;
}
####### include query #######
for($i=1;$i<=$L;$i++)
{
    $am{$it,$i}=substr($queryseq,$i-1,1);
}

####### Henikoff weight $wei{i_seq} ----------->
##### nA{A,i_pos}: number of times A appear at the position:
undef %nA; 
for($i=1;$i<=$Lch;$i++){
    for($j=1;$j<=$it;$j++){
	$nA{$am{$j,$i},$i}++;
    }
}
##### henikoff weight w(i)=sum of 1/rs:
for($i=1;$i<=$it;$i++){
    for($j=1;$j<=$Lch;$j++){
	####### r: number of different residues in j'th position:
	$r=0;
	foreach $A(@AA){
	    $r++ if($nA{$A,$j}>0);
	}
	$A=$am{$i,$j};
	$s2=$nA{$A,$j};
	$w1=1.0/($r*$s2);
	$w{$i}+=$w1;
    }
    $w_all+=$w{$i};
}
#### normalization of w(i):
for($i=1;$i<=$it;$i++){
    $w{$i}/=$w_all;
}
#^^^^^ Henikoff weight finished ^^^^^^^^^^^^^^^^^

########### weighted frequence #################
undef %log;
for($i=1;$i<=$it;$i++){
    for($j=1;$j<=$Lch;$j++){
	$A=$am{$i,$j};
	$log{$j,$A}+=$w{$i};
    }
}
#^^^^^^^^^ Henikoff frequence finished ^^^^^^^^^^^^^

open(freq,">protein.prf");
printf freq "$Lch\n";
for($i=1;$i<=$Lch;$i++){
    printf freq "%3d $seqQ{$i} %3d",$i,$i;
    $norm=0;
    foreach $A(@AA){
        $norm+=$log{$i,$A};
    }
    foreach $A(@AA){
        printf freq "%10.7f",$log{$i,$A}/$norm;
    }
    printf freq "\n";
}
close(freq);

############################### run psipred #########################
print "run psipred...\n";
system("$psipreddir/bin/psipred protein.mtx $psipreddir/data/weights.dat $psipreddir/data/weights.dat2 $psipreddir/data/weights.dat3 >protein.ss");
system("$psipreddir/bin/psipass2 $psipreddir/data/weights_p2.dat 1 1.0 1.0 protein.ss2 protein.ss >protein.horiz");
#`cp protein.ss2 $datadir/protein.ss2`;

$restripletdir="$fudir/lib/ResPre";
system("cd $restripletdir; $python $restripletdir/respre.py $work_dir/deepmsa_protein.aln $work_dir/deepmsa.contact; cd $work_dir");

##### convert it to metapsicov format #######
open(contact,"deepmsa.contact");
open(metacon,">protein.contact");
while($line=<contact>)
{
   @strs=split(/\s/,$line);
   if(abs($strs[0]-$strs[1])>=5)
   {   
   print metacon "$strs[0] $strs[1] 0 8 $strs[2]\n";
   }
}
close(metacon);
close(contact);
#`cp protein.contact $datadir/protein.contact`;
#`cp deepmsa.contact $datadir/respre.contact`;

############### Run CEthreader ###############
system("$bindir/CEdecomposition  -i q -f protein.seq -s protein.horiz -c protein.contact -mtx protein.mtx -o protein.ce -m exp 1.19");
#`cp protein.prf $datadir/protein.prf`;
#`cp protein.mtx $datadir/protein.mtx`;
#`cp protein.ce $datadir/protein.ce`;

print "\n";
print "\n";
print "\n";
print "\n";
print "predicting domain and domain boundary...\n";
############### run FUpred
print "domain boundary is in protein.txt";
#system("$bindir/FUpred -i protein.ce -2c 0.85 -2d 0.66 -chip 15 -label3c 0 -o protein >$work_dir/fu.txt");
#改成输出文件的名字不叫fu.txt,而是叫protein.txt
system("$bindir/FUpred -i protein.ce -2c 0.85 -2d 0.66 -chip 15 -label3c 0 -o protein >$work_dir/protein.txt");
#`cp protein.2c $datadir/`;
#`cp protein.2d $datadir/`;
system("cat $work_dir/fu.txt");
print "\n";
################# endding procedure ######################
$time=`date`;
printf "ending time: $time";
`sync`;
sleep(1);
`rm -fr $work_dir`;
exit();




sub getRandomString
{
    my $len = shift;
    my @W = ('0' .. '9', 'a' .. 'z', 'A' .. 'Z');
    my $str;
    my $i = 0;
    while ($i++ < $len) {
        $str .= $W[rand(@W)];
    }
    return $str;
}

