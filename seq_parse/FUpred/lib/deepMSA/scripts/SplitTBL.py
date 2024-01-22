#!/nfs/amino-library/anaconda/bin/python
import os,sys;
import math;
infilepath=sys.argv[1];
filebasename=os.path.splitext(infilepath)[0];
splitnum=int(sys.argv[2]);
infiler=file(infilepath,'r');
alllines=infiler.readlines();
title=alllines[0:3];
tail=alllines[-10:];
length=math.ceil(1.0*len(alllines[3:-10])/splitnum);
infiler.close()
tmppath='tmp.dat';
tmpfile=file(tmppath,'w');
for allline in alllines[3:-10]:
   tmpfile.write(allline);
tmpfile.close();
infile=file(tmppath,'r');
i=1;
num=1;
line=infile.readline();
templines=line;
while line:
   if(i%length==0):
      print i;      
      outfile=file(filebasename+'.'+str(num)+".tbls",'w');
      for iline in title:
         outfile.write(iline);
      outfile.writelines(templines);
      for jline in tail:
         outfile.write(jline);
      outfile.close();
      templines="";
      num+=1;
   i+=1;
   line=infile.readline();
   templines+=line;

outfile=file(filebasename+'.'+str(num)+".tbls",'w');
for iline in title:
   outfile.write(iline);
outfile.writelines(templines);
for jline in tail:
   outfile.write(jline);
outfile.close();
templines="";
infile.close();
