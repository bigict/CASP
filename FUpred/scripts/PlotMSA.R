args=commandArgs(TRUE);
readmsa<-function(infile)
{
   msafile=file(infile,"r");
   lines=readLines(msafile);
   nrow=length(lines)
   ncol=nchar(lines[1])
   covnum=c();
   for(j in 1:ncol)
   {
      c=0;
      for(i in 1:nrow)
      {
          if(substring(lines[i],j,j)!='-')
          {
              c=c+1;   
          }
      }
      covnum[j]=c;
   }
   close(msafile)
   return(covnum)
}
hhdcovnum=readmsa(args[1])+1
print(hhdcovnum)
#hhpcovnum=readmsa(args[2])+1
jpeg(file=args[2],width=800,height=800)
par(mfrow=c(1,1), mar=c(5,5,1,1))
barplot(hhdcovnum,beside=TRUE,border='grey80',xlab="",ylab="",bty="l",xaxt = "n",yaxt = "n",col='grey80',lwd=5,space=0);
#barplot(hhpcovnum,beside=TRUE,xlab="",ylab="",bty="l",xaxt = "n",yaxt = "n",col='black',lwd=5,ylim=c(1,650),add=T);

mtext("Number of residues per site", side=2, line=3, cex=2.5,font=2);  #### Y
axis(side = 2, lwd = 3,cex.axis=1.5,font=2,mgp=c(0,1.5,0));
mtext("Site",side=1, line=3.5, cex=2.5,font=2);  ##### X
axis(side = 1, lwd = 3,cex.axis=1.5,font=2,mgp=c(0,1.5,0),at=c(1,length(hhdcovnum)),labels=c(1,length(hhdcovnum)));
box(lwd=3,bty="l");

null=dev.off()
