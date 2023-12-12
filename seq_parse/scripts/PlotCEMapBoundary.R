######### Rscript PlotCEMapBoundary.R native.ce predict.ce  out.jpg 1-32;33-50,90-134;51-89; 1-35;36-55,90-134;56-89;
args=commandArgs(TRUE);
library(scales)
readcefile=function(cepath)
{
    ppfile=file(cepath,"r");
    line=readLines(ppfile,n=1);
    tempstr=strsplit(line," ",fixed=TRUE);
    name=tempstr[[1]][4];
    len=as.numeric(tempstr[[1]][5]);
    seq=readLines(ppfile,n=1);
    sec=readLines(ppfile,n=1);
    M=matrix(data=0,nrow=len,ncol=len);
    CSindexs=matrix(data=0,nrow=len,ncol=len);
    contactline=readLines(ppfile,n=1);
    lenc=nchar(contactline);
    contactline=substring(contactline,7,6+lenc);
    contacts=strsplit(contactline,",",fixed=TRUE);
    #print(length(contacts[[1]]))
    for(i in 1:length(contacts[[1]]))
    {
	tempcontactstr=strsplit(contacts[[1]][i]," ",fixed=TRUE);
	indexi=as.numeric(tempcontactstr[[1]][1])+1;
	indexj=as.numeric(tempcontactstr[[1]][2])+1;
        cs=as.numeric(tempcontactstr[[1]][3]);
	M[indexi,indexj]=1;
	M[indexj,indexi]=1;
        CSindexs[indexi,indexj]=cs;
        CSindexs[indexj,indexi]=cs;
    }
    close(ppfile)
    return(list(M,CSindexs,len));
}

PlotBoundary<-function(len,domain_s_vec,domain_e_vec,domain_labels="",mycol="red")
{
    for(i in 1:length(domain_s_vec))
    {   
        d_s=domain_s_vec[i]
        d_e=domain_e_vec[i]
        domain_label=domain_labels[i]
        points(c(d_s,d_s),c(0,len),type='l',lwd=5,col=mycol);
        points(c(d_e,d_e),c(0,len),type='l',lwd=5,col=mycol);
        points(c(0,len),c(d_s,d_s),type='l',lwd=5,col=mycol);
        points(c(0,len),c(d_e,d_e),type='l',lwd=5,col=mycol);
        text(x=(d_s+d_e)/2,y=(d_s+d_e)/2,domain_label,cex=4,col="black",font=2)
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


native=readcefile(args[1])
nativeM=native[[1]]
nativelen=native[[3]]
nxindexs=c();
nyindexs=c();
ncindexs=c();
ncexindexs=c();
for(i in 1:nativelen)
{
	for(j in 1:nativelen)
	{
		nxindexs[length(nxindexs)+1]=i;
		nyindexs[length(nyindexs)+1]=j;
		if(nativeM[i,j]==1)
		{
			ncindexs[length(ncindexs)+1]=alpha("gray80",0.6);
			ncexindexs[length(ncexindexs)+1]=5;
		}
		else
		{
			ncindexs[length(ncindexs)+1]=alpha("white",0.0001);
			ncexindexs[length(ncexindexs)+1]=0;
		}
	}
}


jpeg(args[3],width=2000,height=2000)
par(mfrow=c(1,1), mar=c(15,15,5,5))
plot(nxindexs,nyindexs,col=ncindexs,cex=ncexindexs,pch=20,xlab="",ylab="",bty="o",xaxt = "n",yaxt = "n");


natdomaininfo=ParseBoundary(args[2])
natdomain_s=natdomaininfo[[1]]
natdomain_e=natdomaininfo[[2]]
natdomain_labels=natdomaininfo[[3]]

PlotBoundary(nativelen,natdomain_s,natdomain_e,natdomain_labels,"red")


title(main="",cex.main=8);
mtext("residue index", side=1, line=6, cex=6,font=2)
axis(side = 1, lwd = 5,cex.axis=4,font=2,mgp=c(0,2,0));
mtext("residue index", side=2, line=6, cex=6,font=2)
axis(side = 2, lwd = 5,cex.axis=4,font=2,mgp=c(0,2,0));
box(lwd=5,bty="o");

null=dev.off();
