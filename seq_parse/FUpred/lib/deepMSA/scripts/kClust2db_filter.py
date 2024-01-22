#!/usr/bin/env python
docstring='''
kClust2db.py db.fasta mydb/mydb
    cluster sequences in FASTA file db.fasta using kClust,
    and generate hhblits style database at mydb/mydb

Options:
    -tmpdir=/tmp/$USER/kClust_`date +%N`
        use -tmpdir as temperary folder

    -id=30
        kClust sequence identity cutoff 30%. legal values are: 
        20, 30, 40, 50, 60, 70, 80, 90, 99

    -ncpu=1
        number of CPU threads
'''
import sys,os
import numpy
import subprocess
import shutil
import time
from string import Template

from HHPaths import bin_dict

id2s_dict= { 20:0.52, 30:1.12, 40:1.73, 50:2.33,
    60:2.93, 70:3.53, 80:4.14, 90:4.74, 99:5.28}

kClust_template=Template("mkdir -p $tmpdir/kClust; rm -rf $tmpdir/kClust/*; $kClust -s $s -i $infile -d $tmpdir/kClust -M $m")
kClust_mkAln_template=Template("$kClust_mkAln -c '$clustalo --threads=$ncpu -i $$infile -o $$outfile' -d $tmpdir/kClust --no-pseudo-headers|grep -P '^Filename:'|cut -d' ' -f2")
reformat_template=Template("$reformat fas a3m $filename $tmpdir/a3m/$basename.a3m")
hhblitsdb_template=Template("$hhblitsdb -cpu $ncpu -o $outdb -ia3m $tmpdir/a3m")


def mkdir_if_not_exist(tmpdir):
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

def make_tmpdir(tmpdir):
    ''' create tmp folder '''
    if not tmpdir:
        import random
        tmpdir="/tmp/%s/kClust_%s"%(
            os.getenv("USER"),random.randint(0,10**10))
        while(os.path.isdir(tmpdir)):
            tmpdir="/tmp/%s/kClust_%s"%(
                os.getenv("USER"),random.randint(0,10**10))
    mkdir_if_not_exist(tmpdir)
    sys.stdout.write("created folder %s\n"%tmpdir)
    return tmpdir 

def read_sequences(infile):
    names=[];
    sequences=[];
    fp=file(infile,'r');
    blocks=('\n'+fp.read()).split('\n>');
    for block in blocks:
        lines=block.splitlines();
        if len(lines)>=2:
            names.append(lines[0]);
            sequence='';
            for line in lines[1:]:
                sequence+=line.strip('\n');
            sequences.append(sequence);
    fp.close();
    return names,sequences;


def kClust2db(infile,outdb,tmpdir='.',s=1.12,ncpu=1,m='5000MB'):
    ''' cluster sequences in FASTA file "infile", and generate hhblits
    style database at outdb, return 1'''

    sys.stdout.write("#### cluster input fasta ####\n")
    cmd=kClust_template.substitute(dict(
        kClust=bin_dict["kClust"],
        infile=infile,
        tmpdir=tmpdir,
        m=m,
        s=s,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)

    sys.stdout.write("#### alignment within each cluster ####\n")
    cmd=kClust_mkAln_template.substitute(dict(
        kClust_mkAln=bin_dict["kClust_mkAln"],
        clustalo=bin_dict["clustalo"],
        ncpu=ncpu,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    stdout,stderr=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE).communicate()

    sys.stdout.write("#### reformat fas into a3m ####\n")
    a3mdir=os.path.join(tmpdir,"a3m")
    mkdir_if_not_exist(a3mdir);
    ##### make sure no remains in this folder, or hhblitsdb will generate more than you want!
    os.system("rm "+a3mdir+"/*");
    for filename in stdout.splitlines():
        os.system(reformat_template.substitute(dict(
            reformat=bin_dict["reformat"],
            filename=filename,
            tmpdir=tmpdir,
            basename=os.path.basename(os.path.splitext(filename)[0]),
        )))

    sys.stdout.write("#### build hhblitsdb ####\n")
    mkdir_if_not_exist(os.path.dirname(outdb))
    cmd=hhblitsdb_template.substitute(dict(
        hhblitsdb=bin_dict["hhblitsdb"],
        ncpu=ncpu,
        outdb=outdb,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    return 1;


def kClust2db_pass(infile,outdb,tmpdir='.',s=1.12,ncpu=1,m='5000MB',kclust_timeout=18):
    ''' 
    cluster sequences in FASTA file "infile", and generate hhblits
    style database at outdb, return 1
    if kClust running time is greater than kclust_timeout, then don't build hhblits
    database, return 0.    
    '''

    sys.stdout.write("#### cluster input fasta ####\n")
    k_starttime=time.time();
    cmd=kClust_template.substitute(dict(
        kClust=bin_dict["kClust"],
        infile=infile,
        tmpdir=tmpdir,
        m=m,
        s=s,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)

    k_endtime=time.time();
    k_runningtime=(k_endtime-k_starttime)/60.0;

    if k_runningtime>kclust_timeout:
        return 0;

    sys.stdout.write("#### alignment within each cluster ####\n")
    cmd=kClust_mkAln_template.substitute(dict(
        kClust_mkAln=bin_dict["kClust_mkAln"],
        clustalo=bin_dict["clustalo"],
        ncpu=ncpu,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    stdout,stderr=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE).communicate()

    sys.stdout.write("#### reformat fas into a3m ####\n")
    a3mdir=os.path.join(tmpdir,"a3m")
    mkdir_if_not_exist(a3mdir)
    os.system("rm "+a3mdir+"/*");
    for filename in stdout.splitlines():
        os.system(reformat_template.substitute(dict(
            reformat=bin_dict["reformat"],
            filename=filename,
            tmpdir=tmpdir,
            basename=os.path.basename(os.path.splitext(filename)[0]),
        )))

    sys.stdout.write("#### build hhblitsdb ####\n")
    mkdir_if_not_exist(os.path.dirname(outdb))
    cmd=hhblitsdb_template.substitute(dict(
        hhblitsdb=bin_dict["hhblitsdb"],
        ncpu=ncpu,
        outdb=outdb,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    return 1;



def filter_kclust(list_in,list_out,Nmax):
    '''
    list_in      - path for kclust and kclust_mkaln output list 
    list_out     - path for filtered clusters path list
    Nmax         - Max numbers of clusters
    '''
    
    infile=file(list_in,'r');
    
    lines=infile.readlines();

    clusters_path_list=[];
    cscores=[];
    for line in lines:
        cluster_path=line.strip('\n');
        (names,seqs)=read_sequences(cluster_path);
        
        cluster_score='';
        evs=[];
        for name in names:
            strs=name.split('#');
            for tmpstr in strs:
                if tmpstr.startswith("psi-ev:"):
                    e=tmpstr.split(':')[1];
                    if e.startswith('e'):
                        e='1'+e;
                    evs.append(float(e));
        if len(evs)==0:
            cluster_score=0;
        else:
            cluster_score=min(evs);   ##### need check #####
        
        cscores.append(cluster_score);
        clusters_path_list.append(cluster_path);        
    infile.close();    
    ############# sort clusters and take top Nmax clusters ###############
    filtered_paths=[];
    cslist=numpy.array(cscores);
    indexlist = numpy.argsort(cslist);
    num=0;
    for index in indexlist:
        filtered_paths.append(clusters_path_list[index]);
        num+=1;
        if num>=Nmax:
            break;

    outfile=file(list_out,'w');
    for filtered_path in filtered_paths:
         outfile.write(filtered_path+'\n');
    outfile.close();
    print "filter cluster number from %d to %d"%(len(clusters_path_list),len(filtered_paths));
    return
    

def kClust2db_filter(infile,outdb,Nmax=999999,tmpdir='.',s=1.12,ncpu=1,m='5000MB'):
    ''' cluster sequences in FASTA file "infile", and generate hhblits
    style database at outdb
    the infile should have psiblast e-value, then basing on e-value 
    filter cluster num
    '''

    sys.stdout.write("#### cluster input fasta ####\n")
    
    cmd=kClust_template.substitute(dict(
        kClust=bin_dict["kClust"],
        infile=infile,
        tmpdir=tmpdir,
        m=m,
        s=s,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    
    
    sys.stdout.write("#### alignment within each cluster ####\n")
    cmd=kClust_mkAln_template.substitute(dict(
        kClust_mkAln=bin_dict["kClust_mkAln"],
        clustalo=bin_dict["clustalo"],
        ncpu=ncpu,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    stdout,stderr=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE).communicate()
    
    allpath=tmpdir+'/kcluster.path';
    
    allfile=file(allpath,'w');
    allfile.write(stdout);
    allfile.close();
    
    filterpath=tmpdir+'/kcluster_filter.path';

    filter_kclust(list_in=allpath,list_out=filterpath,Nmax=Nmax);
    filterfile=file(filterpath,'r');
    filtercontext=filterfile.read();
    filterfile.close();

    sys.stdout.write("#### reformat fas into a3m ####\n")
    a3mdir=os.path.join(tmpdir,"a3m")
    mkdir_if_not_exist(a3mdir)
    os.system("rm "+a3mdir+"/*");
    for filename in filtercontext.splitlines():
        os.system(reformat_template.substitute(dict(
            reformat=bin_dict["reformat"],
            filename=filename,
            tmpdir=tmpdir,
            basename=os.path.basename(os.path.splitext(filename)[0]),
        )))

    sys.stdout.write("#### build hhblitsdb ####\n")
    mkdir_if_not_exist(os.path.dirname(outdb))
    cmd=hhblitsdb_template.substitute(dict(
        hhblitsdb=bin_dict["hhblitsdb"],
        ncpu=ncpu,
        outdb=outdb,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    return 1;

def kClust2db_filter_pass(infile,outdb,Nmax=999999,tmpdir='.',s=1.12,ncpu=1,mem='5000MB',kclust_timeout=18):
    ''' 
    cluster sequences in FASTA file "infile", and generate hhblits
    style database at outdb, return 1.
    the infile should have psiblast e-value, then basing on e-value 
    filter cluster num, if kClust running time is greater than kclust_timeout, then don't build hhblits 
    database, return 0.
    '''

    sys.stdout.write("#### cluster input fasta ####\n")
    k_starttime=time.time();
    cmd=kClust_template.substitute(dict(
        kClust=bin_dict["kClust"],
        infile=infile,
        tmpdir=tmpdir,
        m=mem,
        s=s,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    k_endtime=time.time();
    k_runningtime=(k_endtime-k_starttime)/60.0;

    if k_runningtime>kclust_timeout:
        return 0;
    
    sys.stdout.write("#### alignment within each cluster ####\n")
    cmd=kClust_mkAln_template.substitute(dict(
        kClust_mkAln=bin_dict["kClust_mkAln"],
        clustalo=bin_dict["clustalo"],
        ncpu=ncpu,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    stdout,stderr=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE).communicate()
    
    allpath=tmpdir+'/kcluster.path';
    
    allfile=file(allpath,'w');
    allfile.write(stdout);
    allfile.close();
    
    filterpath=tmpdir+'/kcluster_filter.path';

    filter_kclust(list_in=allpath,list_out=filterpath,Nmax=Nmax);
    filterfile=file(filterpath,'r');
    filtercontext=filterfile.read();
    filterfile.close();

    sys.stdout.write("#### reformat fas into a3m ####\n")
    a3mdir=os.path.join(tmpdir,"a3m")
    mkdir_if_not_exist(a3mdir)
    os.system("rm "+a3mdir+"/*");
    for filename in filtercontext.splitlines():
        os.system(reformat_template.substitute(dict(
            reformat=bin_dict["reformat"],
            filename=filename,
            tmpdir=tmpdir,
            basename=os.path.basename(os.path.splitext(filename)[0]),
        )))

    sys.stdout.write("#### build hhblitsdb ####\n")
    mkdir_if_not_exist(os.path.dirname(outdb))
    cmd=hhblitsdb_template.substitute(dict(
        hhblitsdb=bin_dict["hhblitsdb"],
        ncpu=ncpu,
        outdb=outdb,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    return 1;




if __name__=="__main__":
    seqID=30
    ncpu=1
    tmpdir=''
    
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-id="):
            seqID=float(arg[len("-id="):])
            if seqID<1:
                seqID=100*seqID
            seqID=int(seqID)
        elif arg.startswith("-ncpu="):
            ncpu=int(arg[len("-ncpu="):])
        elif arg.startswith("-tmpdir="):
            tmpdir=os.path.abspath(arg[len("-tmpdir="):])
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! No such option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if not seqID in id2s_dict:
        sys.stderr.write("ERROR! Illegal sequence identity cutoff %d\n"%seqID)
        exit()
    s=id2s_dict[seqID]

    if len(argv)!=2:
        sys.stderr.write(docstring)
        exit()

    infile=os.path.abspath(argv[0])
    outdb=os.path.abspath(argv[1])
    tmpdir=make_tmpdir(tmpdir)

    kClust2db_filter(infile,outdb,300,tmpdir,s,ncpu)
    #if os.path.isdir(tmpdir):
    #    shutil.rmtree(tmpdir)
