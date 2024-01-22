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

def refomarta3m(a3mfile,outbase):
    
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
    Nmsa=len(blocks[k:]);
    #print queryblock
    queryblock='>'+queryblock;
    outfile=file(outbase,'w');
    querylines=queryblock.split('\n');
    strline='';
    for queryline in querylines[1:]:
        strline+=queryline;
    querylen=len(strline);
    outfile.write("%5d %5d protein\n"%(Nmsa,querylen));
    for queryline in querylines[1:]:
        if len(queryline)==50:
            outfile.write(queryline+'*\n');
        elif len(queryline)>0:
            outfile.write(queryline+'\n');

    
    homoblocks=blocks[k+1:]
    if querylen>1:
        for homoblock in homoblocks:
            blockstrs=homoblock.split('\n');
        
            for line in blockstrs[1:]:
                if len(line)>0:
                    outfile.write(line+'\n')
        
    
    a3m.close();
    outfile.close();
    return len(homoblocks);

if __name__=="__main__":

    infile=sys.argv[1];    ####a3m file
    outbase=sys.argv[2];

    refomarta3m(infile,outbase);
