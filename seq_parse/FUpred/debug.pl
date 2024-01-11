#!/usr/bin/perl
use Math::Trig;
use File::Basename;
use Cwd 'abs_path';

######## set your variables here
#$msa_hhblitsdb="/nfs/amino-home/zhng/local_library/uniclust30_2017_04/uniclust30_2017_04";
$msa_hhblitsdb="/public/dbu/data/af2_db/uniclust30/uniclust30_2018_08/uniclust30_2018_08";
#$msa_jackhmmerdb="/nfs/amino-home/zhng/local_library/uniref90_04042018/uniref90.fasta";
$msa_jackhmmerdb="/public/dbu/data/fupred_db/uniref90_db/uniref90.fasta";
#$msa_hmmsearchdb="/nfs/amino-home/zhng/local_library/metaclust/metaclust.fasta";
$msa_hmmsearchdb="/public/dbu/data/fupred_db/metaclust_db/metaclust_2017_05.fasta";
#$python="/nfs/amino-home/zhng/local_library/anaconda3/bin/python";
$python="/public/dbu/local/anaconda3/envs/respre/bin/python";
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


########### get mtx and prf from deepmsa ######################

system("$python2 $hhsuite/scripts/build_MSA.py -ncpu=$cpunum protein.seq -hhblitsdb=$msa_hhblitsdb -jackhmmerdb=$msa_jackhmmerdb  -hmmsearchdb=$msa_hmmsearchdb -tmpdir=$work_dir/deepMSA");
#用echo $?查看上一条命令的返回值，如果返回值为0，则表示执行成功，否则执行失败，打印出错误信息。
if($?){
    print "build_MSA.py failed\n";
    exit();
}
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


################# endding procedure ######################
$time=`date`;
printf "ending time: $time";
`sync`;
sleep(1);
#`rm -fr $work_dir`;
exit();




sub getRandomString
{
    # my $len = shift;
    # my @W = ('0' .. '9', 'a' .. 'z', 'A' .. 'Z');
    # my $str;
    # my $i = 0;
    # while ($i++ < $len) {
    #     $str .= $W[rand(@W)];
    # }
    #令str值为rxI2n5KBC6
    my $str="rxI2n5KBC6";
    return $str;
}

