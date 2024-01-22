#!/bin/env python
import os;
import sys;
import commands;
import HHPaths;
import numpy
import shutil
import subprocess
import time
from string import Template,ascii_lowercase
from HHPaths import bin_dict

hhlib=os.environ["HHLIB"]
# format db for psiblast
# $infile      - input fasta file
# $mydbtitle   - formatted database name
formatdb_template=Template(bin_dict["formatdb"]+" -n $mydbtitle -i $infile -p T")

# run psiblast 
# $infile      -  input query sequence
# $ncpu        -  number of cpu cores
# $noutput     -  numbers of output
# $evalue      -  e-value cut off
# $db          -  database
# $blastout    -  psiblast output file
# $iter        -  psiblast iteration times
psiblast_template=Template(bin_dict["psiblast"]+" -a $ncpu -b $noutput -m $m -j $iter -h $evalue -d $db -i $infile -o $blastout")


def a3m2querydb(a3mfile,outbase):
    
    a3m=file(a3mfile,'r');
    ###### split to blocks #####
    content=a3m.read();
    content='\n'+content;
    blocks=content.split('\n>');
    
    #query blocks[1]
    for k in range(1,len(blocks)):
        if not blocks[k].startswith('ss_'):
            break;
    queryblock=blocks[k];
    #print queryblock
    queryblock='>'+queryblock;
    queryfile=file(outbase+'-query.fasta','w');
    queryfile.write(queryblock+'\n')
    queryfile.close();
    
    dbfile=file(outbase+'-db.fasta','w');
    homoblocks=blocks[k+1:]
    for homoblock in homoblocks:
        blockstrs=homoblock.split('\n');
        id='>'+blockstrs[0];
        dbfile.write(id+'\n');
        for line in blockstrs[1:]:
            line=line.replace('-','');
            line=line.replace('.','');
            line=line.upper();
            dbfile.write(line)
        dbfile.write('\n')
    dbfile.close();
    a3m.close();
    return len(homoblocks);

def mkblastout(queryfile,dbfasta,noutput,outbase):
    ####### formatdb
    cmd=formatdb_template.substitute(
        infile     =dbfasta,
        mydbtitle  ="mydb",
    )
    sys.stdout.write(cmd+'\n')
    os.system(cmd)

    m=0;
    iter=3;
    ncpu=1;
    evalue=0.005;
    db="mydb";
    blastout=outbase+'-blast.out'
    '''  run psiblast with evalue[default 0.001] '''
    cmd=psiblast_template.substitute(
        ncpu       =ncpu,
        noutput    =noutput,
        iter       =iter,
        evalue     =evalue,
        db         =db,
        infile     =queryfile,
        blastout   =blastout,
        m          =m,
    )
    sys.stdout.write(cmd+'\n')
    os.system(cmd)



if __name__=="__main__":

    infile=sys.argv[1];    ####a3m file
    outbase=sys.argv[2];

    ####### reformat a3m
    nout=a3m2querydb(infile,outbase);
    
    ##### output outbase-query.fasta and outbase-db.fasta
    mkblastout(outbase+"-query.fasta",outbase+"-db.fasta",nout,outbase);

    ##### remove files
    #os.system('rm '+outbase+'-query.fasta');
    os.system('rm '+outbase+'-db.fasta');    
    os.system('rm mydb.*')

