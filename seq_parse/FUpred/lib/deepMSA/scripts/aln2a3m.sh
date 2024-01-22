#!/bin/sh
alnfile=$1
a3mfile=$2
i=0;
while read line
do
if [ "$i" -eq "0" ]
then
   echo ">protein" >$a3mfile
   echo "$line"    >>$a3mfile
   i=`expr $i + 1`
else
   echo ">seq_$i"  >>$a3mfile
   echo "$line"   >>$a3mfile
   i=`expr $i + 1`
fi
done<$alnfile
