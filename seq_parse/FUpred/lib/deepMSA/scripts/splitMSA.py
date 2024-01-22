#!/bin/env python
### example use: python splitfasta.py ex/deepmsa.fas ex/subfiles/ 10
###
###
### input: 1. fastafile: alignments of query and multiple templates
###        2. outdir: place the output files
###        3. number of splitted files
###
###
### output: 1. splitted subfastafiles (output to outdir): alignments of query and multiple templates
###
###
### description: the sequence length of input fastafile is long, split the input fastafile to several subfastafiles. The sequence length of each subfastafile = The total sequence length of input fastafile / numbe of splitted files. If the total sequence length of input fastafiel cannot be divided with no remainder, the last outputfile will be more or less but other outputfiles are exactly the same number.
###              Eg. The input fastafile has a length of 303. We want to split it into 3 subfastafiles. So the number of sequence length for each subfastafile is 303/3=101. Finally, alignments of 1-101 is the first subfastafile, alignments of  102-202 is the second subfastafile, and 203-303 is the third subfastafile. 


import sys;
import os;

def splitfasta(infastafile,outdir,nfiles):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # read input fastafile to a list
    fastaseqlist = [] # store the "order" of fasta seqs in input fastafile
    fastanamelist = [] # store the "order" of fasta names in input fastafile
    pointer = "" # store fasta seq

    ifr = open(infastafile,'r')
    for line in ifr:
        if line.strip().startswith('>'):
            if pointer != "":
                fastaseqlist.append(pointer)
            pointer = ""
            fastanamelist.append(line.strip()[1:])
        else:
            pointer = pointer + line.strip()
    fastaseqlist.append(pointer)
    ifr.close()

    # sequence length
    seqlen = len(fastaseqlist[0])
    subseqlen = seqlen/nfiles

    # write outputfile 
    newfastanamelist = [fastanamelist[0]] +  range(1, len(fastanamelist))
    for j in range(1, nfiles):
        outfile = outdir + '/sub_' + str(j) + '.fas'
        ofw = open(outfile, 'w')
        for i in range(0, len(fastanamelist)):
            ofw.write(">" + str(newfastanamelist[i]) + "#" + str((j-1)*subseqlen + 1) + "-" + str(j*subseqlen) + "\n")
            ofw.write(fastaseqlist[i][(j-1)*subseqlen:j*subseqlen] + "\n")    
        ofw.close()

    # write the last outputfile
    outfile = outdir + '/sub_' + str(nfiles) + '.fas'
    ofw = open(outfile, 'w')
    for i in range(0, len(fastanamelist)):
        ofw.write(">" + str(newfastanamelist[i]) + "#" + str((nfiles-1)*subseqlen + 1) + "-" + str(seqlen) + "\n")
        ofw.write(fastaseqlist[i][(nfiles-1)*subseqlen:] + "\n")
    ofw.close()


if  __name__ == "__main__":
    infastafile = sys.argv[1]
    outdir = sys.argv[2]
    nfiles = int(sys.argv[3])
    splitfasta(infastafile,outdir,nfiles)
