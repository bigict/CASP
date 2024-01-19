                                                         FUpred version 1.0
                                                       ============================

Directories in this source folder:

    bin (CEdecomposition  and FUpred)
    lib (include dssp, ResPRE, deepMSA and psipred4)
    scripts (PlotCEMapBoundary.R,  PlotMSA.R and PlotScore.R are visualization tools)
    example (example)
    run-FUpred.pl (main program)
    readme.txt
    licence.txt


Overview

FUpred is a contact map-based domain prediction method which utilizes a recursion strategy 
to detect domain boundary based on predicted contact-map and secondary structure information. 
Large scale benchmark analysis shows that FUpred has significantly better ability of domain
boundary prediction than threading-based method and machine learning-based methods. 
Particularly, our method has obviously excellent performance in detecting discontinuous domain 
boundary than current methods.

Reference

FUpred: Detecting protein domains through deep-learning based contact map prediction. Bioinformatics, (2020)
W Zheng, X.G. Zhou, Q.Q.G. Wuyun, R Pearce, Y Li and Y Zhang.

###############################################################################################################################################
1. Installation

   You can directly run it in Linux system, but you need change the variables in run-FUpred.pl following the instruction below.


################################################################################################################################################
2. Programs instruction

The main program is run-FUpred.pl, given a fasta format protein sequence file
you can run  

./run-FUpred.pl sequence.fasta

where sequence.fasta is your input file.
you can try

./run-FUpred.pl ./example/6paxa.fasta

the predicted results will be output to screen like this:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
predicting domain and domain boundary...
domain boundary is:1-65;66-133;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if you want to customize the FUpred program with different parameters, then you can read (2.1)
and (2.2) for help, else you can skip them.

-------------------------------------------------------------------------------------------------------------------
(2.1). CEdecomposition

This program is contact map eigen decomposition method, design for FUpred/CEthreading program.

Usage: CEdecomposition -i inputtype[q: query t: template] -f fasta [-n native pdb] -s psipred.horiz [-d dssp] -c metapsicov_contact_format[if -n no predicted_contact need] -o outfile -m [linear num or exp num or top num] -mtx psiblastmtx

opitions:
          -i: input type q: query                  [design for CEthreading]
                         t: template               [design for CEthreading]
                         qnm: query without mtx    [design for FUpred]
                         tnm: template without mtx [design for FUpred]
for query (q and qnm):
          -f: input fasta for query
          -s: psipred horiz out file
          -c: metapsicov format contact file
          -m choose one of following three cutoff, default exp
          linear: linear model for top num*L contact cutoff, default num=2
          exp: exponent model for top L^num contact cutoff, default num=1.2
          top: top num[fixed] contact cutoff
for template (t and tnm):
          -n: native structure for template
          -d: dssp file
q and t common:
          -mtx: psiblast mtx file
          -o: output file
qnm and tnm common:
          -o: output file


How to build ce file by you own seleced contacts

For example:
original contact map by CASP format

1 19 0 8 0.991
1 18 0 8 0.71
91 103 0 8 0.700
...........
32 54 0 8 0.001

(total 3000 contacts)

You can select any contacts as you want, for example (only two contacts and ignore the confidence scores)

91 103 0 8 0.700
32 54 0 8 0.001

Then write these to a file (mycontact.con)
Then use CEdecomposition do eigen decomposition

CEdecomposition -i qnm  -f fastafile  -s psipred.horiz  -c mycontact.con -o outfile -m  top 2 

Then you will get a input file basing on a conatct map only contains two contacts. 
This "-m top" parameters are useful when you build your own contact map. You don't need change any source code!

-------------------------------------------------------------------------------------------------------------------
(2.2). FUpred

This is the contact map-based recursion strategy domain partition program which uses ce format file as input.

Usage:  -i inputfile [xxxx.ce format ] -2c [two continuous domain cutoff] -2d [two discontinuous domain cutoff] -chip [chip length] -label3c [use 3c or not]

explanation:
       -chip 
                     when split the squence to domains using recursion strategy, the protein will be 
                     split into small fragments, if the length of the fragment is less than this chip
                     we will merge the fragment to last stage domain fragment to avoid too many small 
                     fragments in final results.
       -label3c
                     this is the parameter that indicate whether FUpred will use a 3-domain continuous
                     domain partition scoring function to refine the domain bounddary. 0 means not use it
                     1 means use it.  

##################################################################################################################################
3. Databases:

(1). For build deep MSA

User must download three sequence databases

Unicluset30 from http://gwdu111.gwdg.de/~compbiol/uniclust/2018_08/ (notice you must use hhsuite-2.0 format databases, not fasta format)
Uniref90    from https://www.uniprot.org/downloads  (fasta format)
Metaclust   from https://metaclust.mmseqs.org/2018_06/  (fasta format)

To use DeepMSA, you need download uniclust30, uniref90 and metaclust from 
http://gwdu111.gwdg.de/~compbiol/uniclust/2017_04/uniclust30_2017_04_hhsuite.tar.gz , 
ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz , and 
https://metaclust.mmseqs.org/2017_05/metaclust_2017_05.fasta.gz. 
After you unpack them, then use $pkgdir/lib/deepMSA/bin/esl-sfetch to create .ssi index 
for uniref90 and metaclust, here $pkgdir means the path where you put the FUpred package. 
For example, if the uniref90 database in uniref90 folder is named as uniref90.fasta, 
then go to uniref90 folder, run $pkgdir/lib/deepMSA/bin/esl-sfetch --index uniref90.fasta, 
you will find a new file named as uniref90.fasta.ssi after the command done. Then do the same thing to metaclust database.

And change variables "$msa_hhblitsdb" "$msa_jackhmmerdb" "$msa_hmssearchdb" in run-FUpred.pl
 
########################################################################################################################################################
4. Dependence

We attach "psipred-4.01", "DeepMSA" and "ResPRE" in lib folder, to use ResPRE, please read README file 
in ResPRE folder and install python with "numpy", "scipy" and "pytorch". then change variable "$python" in run-FUpred.pl
You also need change the $python2 to a python with version 2.7 (suggest Anaconda2).
Here DeepMSA requires python2 and ResPRE requires python3, we recommands you install Anaconda2 and Anaconda3.
#########################################################################################################################################################


If you have any question, please contact with Wei Zheng (zhengwei@umich.edu or jlspzw139@sina.com).
