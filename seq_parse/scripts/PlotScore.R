library('scales')
library('plot3D')
args=commandArgs(TRUE);

PlotBoundary2d<-function(len,domain_s_vec,domain_e_vec,domain_labels="",mycol="red")
{
    for(i in 1:length(domain_s_vec))
    {   
        d_s=domain_s_vec[i]
        d_e=domain_e_vec[i]
        domain_label=domain_labels[i]
        points(c(d_s,d_s),c(0,len),type='l',lwd=5,col=mycol,lty=3);
        points(c(d_e,d_e),c(0,len),type='l',lwd=5,col=mycol,lty=3);
        points(c(0,len),c(d_s,d_s),type='l',lwd=5,col=mycol,lty=3);
        points(c(0,len),c(d_e,d_e),type='l',lwd=5,col=mycol,lty=3);
        #text(x=(d_s+d_e)/2,y=(d_s+d_e)/2,domain_label,cex=4,col="black",font=2)
    }
}
PlotBoundary2c<-function(len,domain_s_vec,domain_e_vec,domain_labels="",mycol="red")
{
    for(i in 1:length(domain_s_vec))
    {   
        d_s=domain_s_vec[i]
        d_e=domain_e_vec[i]
        domain_label=domain_labels[i]
        points(c(d_s,d_s),c(0,len),type='l',lwd=5,col=mycol,lty=3);
        points(c(d_e,d_e),c(0,len),type='l',lwd=5,col=mycol,lty=3);
        #points(c(0,len),c(d_s,d_s),type='l',lwd=5,col=mycol);
        #points(c(0,len),c(d_e,d_e),type='l',lwd=5,col=mycol);
        #text(x=(d_s+d_e)/2,y=(d_s+d_e)/2,domain_label,cex=4,col="black",font=2)
    }
}

ParseBoundary<-function(domainstring)
{
    domain_s=c()
    domain_e=c()
    domain_labels=c()
    domainstrs=strsplit(domainstring,";",fixed=TRUE)
    maxd=0;
    j=0
    for(i in 1:length(domainstrs[[1]]))
    {
         if(domainstrs[[1]][i]!="")
         {
            maxd=maxd+1
            if(grepl(",",domainstrs[[1]][i],fixed=TRUE))
            {
                subdomains=strsplit(domainstrs[[1]][i],",",fixed=TRUE)
                for(k in 1:length(subdomains[[1]]))
                {
                    j=j+1;
                    domain_labels[j]=paste(paste("D",maxd,sep=""),k,sep="-");
       	       	    domaintmp=strsplit(subdomains[[1]][k],"-",fixed=TRUE)
       	       	    domain_s[j]=as.numeric(domaintmp[[1]][1])
                    domain_e[j]=as.numeric(domaintmp[[1]][2])

                }
            }
            else
            {
                j=j+1
                domain_labels[j]=paste("D",maxd,sep="");
                domaintmp=strsplit(domainstrs[[1]][i],"-",fixed=TRUE)
                domain_s[j]=as.numeric(domaintmp[[1]][1])
       	       	domain_e[j]=as.numeric(domaintmp[[1]][2])
            }
         }
    }
    return(list(domain_s,domain_e,domain_labels))
}


DCscores2c=c()

score2cfile=read.table(args[1],header=F,sep=" ")
DCscores2c=score2cfile$V2;
psedo=4;
for(k in 1:length(DCscores2c))
{
    if(is.nan(DCscores2c[k])||is.infinite(DCscores2c[k]))
    {
       DCscores2c[k]=psedo;
    }
	else
	{
		if(DCscores2c[k]>psedo)
		{
			DCscores2c[k]=psedo;
		}
	}
}
len=length(DCscores2c)
DCscores2c[length(DCscores2c)]=psedo;
score2dfile=read.table(args[2],header=F,sep=" ")
DCscores2d=matrix(data=psedo,nrow=len,ncol=len);
for(i in 1:len)
{
   for(j in 1:len)
   {
		DCscores2d[i,j]=as.numeric(score2dfile[i,j])
   }
}
#DCscores2c

for(i in 1:length(DCscores2c))
{
	for(j in 1:length(DCscores2c))
	{
		if(is.nan(DCscores2d[i,j])||is.infinite(DCscores2d[i,j]))
        {
            DCscores2d[i,j]=psedo;
        }
		else
		{
			if(DCscores2d[i,j]>psedo)
			{
				DCscores2d[i,j]=psedo;
			}
		}
	}
}
#DCscores2d
natdomaininfo=ParseBoundary(args[3])
natdomain_s=natdomaininfo[[1]]
natdomain_e=natdomaininfo[[2]]
natdomain_labels=natdomaininfo[[3]]



jpeg(args[4],width=4000,height=2000)
par(mfrow=c(1,2), mar=c(10,12,5,5))

plot(seq(1,len,1),DCscores2c,col='black',type='l',lwd=5,xlim=c(0,len),ylim=c(0,psedo),pch=20,xlab="",ylab="",bty="l",xaxt = "n",yaxt = "n");

mtext("Residue index", side=1, line=7, cex=5,font=2)
axis(side = 1, lwd = 5,cex.axis=4,font=2,mgp=c(0,3,0));
mtext("FU-score", side=2, line=7, cex=5,font=2)
axis(side = 2, lwd = 5,cex.axis=4,font=2,mgp=c(0,3,0));
PlotBoundary2c(len,natdomain_s,natdomain_e,natdomain_labels,"red")
box(lwd=5,bty="l");

image2D(x=seq(1,len,1),y=seq(1,len,1),DCscores2d,shade = 0.2, rasterImage = TRUE,contour = list(col = "black", labcex = 6, lwd = 3,dist=10, alpha = 0.5),alpha=0.75,colkey = FALSE,xlab="",ylab="",xaxt = "n",yaxt = "n")
mtext("Residue index", side=1, line=7, cex=5,font=2)
axis(side = 1, lwd = 5,cex.axis=4,font=2,mgp=c(0,3,0));
mtext("Residue index", side=2, line=7, cex=5,font=2)
axis(side = 2, lwd = 5,cex.axis=4,font=2,mgp=c(0,3,0));
PlotBoundary2d(len,natdomain_s,natdomain_e,natdomain_labels,"black")
box(lwd=5,bty="l");

null=dev.off()