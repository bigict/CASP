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
print "$infile\n";
#$datadir=$ENV{'PWD'};
#$datadir=abs_path($ARGV[0]);
#$infile="$datadir/seq.fasta";
$scriptsdir=dirname(abs_path(__FILE__));
$fupredbindir="$scriptsdir/bin";
$bindir=abs_path("$scriptsdir/./bin");
$fudir=abs_path("$bindir/../");
#print "$bindir\n";
$method="respre";


my $random = &getRandomString(10);
$work_dir = "/tmp/$ENV{'USER'}/FUD_$random";
$deepmsa_tmp_dir = "/tmp/$ENV{'USER'}/DeepMSAtmp_$random";
$inputdir = "/tmp/$ENV{'USER'}/Input_$random";

# $work_dir="/tmp/$ENV{'USER'}/FUD_".&getRandomString(10);
# $deepmsa_tmp_dir="/tmp/$ENV{'USER'}/DeepMSAtmp_".&getRandomString(10);
# $inputdir="/tmp/$ENV{'USER'}/Input_".&getRandomString(10);
print "work_dir is:$work_dir\n";
print "inputdir is:$inputdir\n";
print "deepmsa_tmp_dir is:$deepmsa_tmp_dir\n";
$psipreddir="$fudir/lib/psipred-4.01";
$ENV{'HHLIB'}="$fudir/lib/deepMSA";
$hhsuite=$ENV{'HHLIB'};
print "$hhsuite\n";
################ working directory ########################
`/bin/mkdir -p $work_dir`;
chdir "$work_dir";
`/bin/rm -rf $work_dir/*`;
`/bin/mkdir -p $inputdir`;
`/bin/rm -rf $inputdir/*`;
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
open(seq,">","$inputdir/protein.seq") or die "can't open $inputdir/protein.seq\n";
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

printf "start running deepmsa\n";

system("$python2 $hhsuite/scripts/build_MSA.py -ncpu=$cpunum $inputdir/protein.seq -hhblitsdb=$msa_hhblitsdb -jackhmmerdb=$msa_jackhmmerdb  -hmmsearchdb=$msa_hmmsearchdb -tmpdir=$work_dir/deepMSA -outdir=$deepmsa_tmp_dir");
#用echo $?查看上一条命令的返回值，如果返回值为0，则表示执行成功，否则执行失败，打印出错误信息。
my $deepmsa_exit_code = $? >> 8;
print "deepmsa exit code:$deepmsa_exit_code\n";


my $current_directory = abs_path();
print "Current Directory: $current_directory\n";
`cp $deepmsa_tmp_dir/protein.aln $deepmsa_tmp_dir/deepmsa_protein.aln`;
`cp $deepmsa_tmp_dir/protein.a3m $deepmsa_tmp_dir/deepmsa_protein.a3m`;
#check whether cp is successful
my $cp_exit_code = $? >> 8;
print "cp_exit_code: $cp_exit_code\n";



`$ENV{'HHLIB'}/scripts/a3m2psiblast_v2.pl $deepmsa_tmp_dir/deepmsa_protein.a3m psitmp -a3m -neff 7 -output $deepmsa_tmp_dir`;   ## for ss

`cp $deepmsa_tmp_dir/psitmp.mtx  $deepmsa_tmp_dir/protein.mtx`;
`cp $deepmsa_tmp_dir/psitmp.pssm $deepmsa_tmp_dir/protein.pssm`;
my $cp_exit_code = $? >> 8;
print "cp_exit_code of psitmp: $cp_exit_code\n";


open(msa,"$deepmsa_tmp_dir/deepmsa_protein.aln");
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

open(freq,">","$inputdir/protein.prf") or die "can't open $inputdir/protein.prf\n";
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
system("$psipreddir/bin/psipred $deepmsa_tmp_dir/protein.mtx $psipreddir/data/weights.dat $psipreddir/data/weights.dat2 $psipreddir/data/weights.dat3 > $inputdir/protein.ss");
my $psipred_exit_code = $? >> 8;
print "psipred exit code:$psipred_exit_code\n";


system("$psipreddir/bin/psipass2 $psipreddir/data/weights_p2.dat 1 1.0 1.0 protein.ss2 $inputdir/protein.ss > $inputdir/protein.horiz");
my $psipass2_exit_code = $? >> 8;
print "psipass2 exit code:$psipass2_exit_code\n";

#`cp protein.ss2 $datadir/protein.ss2`;
print "===========start run respre=============\n";
$restripletdir="$fudir/lib/ResPre";
system("cd $restripletdir; $python $restripletdir/respre.py $deepmsa_tmp_dir/deepmsa_protein.aln $work_dir/deepmsa.contact; cd $work_dir");
my $restripletdir_exit_code = $? >> 8;
print "restripletdir exit code:$restripletdir_exit_code\n";

##### convert it to metapsicov format #######
open(contact,"deepmsa.contact") or die "can't open deepmsa.contact\n";
open(metacon,">protein.contact") or die "can't open protein.contact\n";
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
# system("cd $scriptsdir");
# my cd_exit_code=$? >> 8;
# print "cd_exit_code: $cd_exit_code\n";
`pwd`;

#print "scriptsdir:$scriptsdir\n";

system("ls -lht");

system("$bindir/CEdecomposition  -i q -f $inputdir/protein.seq -s $inputdir/protein.horiz -c $work_dir/protein.contact -mtx $deepmsa_tmp_dir/protein.mtx -o $inputdir/protein.ce -m exp 1.19");
my $CEdecomposition_exit_code = $? >> 8;
print "CEdecomposition exit code:$CEdecomposition_exit_code\n";

#`cp protein.prf $datadir/protein.prf`;
#`cp protein.mtx $datadir/protein.mtx`;
#`cp protein.ce $datadir/protein.ce`;

print "\n";
print "\n";
print "\n";
print "\n";
print "predicting domain and domain boundary...\n";
############### run FUpred
my($protein_name,$protein_path,$protein_suffix)=fileparse($infile,qr/\.[^.]*/);
print("protein_name: $protein_name\n");
$output_dir="$scriptsdir/../FUpredOutput";
print("output_dir: $output_dir\n");
`/bin/mkdir -p $output_dir`;


system("$bindir/FUpred -i $inputdir/protein.ce -2c 0.85 -2d 0.66 -chip 15 -label3c 0 -o protein >$work_dir/fu.txt");
#system("$bindir/FUpred -i protein.ce -2c 0.85 -2d 0.66 -chip 15 -label3c 0 -o protein >$work_dir/protein-FUpred.txt");
#`cp protein.2c $datadir/`;
#`cp protein.2d $datadir/`;
print("domain boundary is:");
system("cat $work_dir/fu.txt");
print "\n";
print "domain boundary is in $output_dir/$protein_name-FUpred.txt  \n";
system("/bin/mv $work_dir/fu.txt $output_dir/$protein_name-FUpred.txt");
################# endding procedure ######################
$time=`date`;
printf "ending time: $time";
`sync`;
sleep(1);
#`rm -fr $work_dir`;
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

