library(pegas)
library(haplotypes) 
path = "~/mafft_alignment"
setwd(path)
cladeArray<-c("A","B","C","D","E","F","G","H","I")
#cladeArray，Loop processing cladeA-I：
#cladeA.mafft.aln.fas
#cladeB.mafft.aln.fas
#cladeC.mafft.aln.fas
#cladeD.mafft.aln.fas
#cladeE.mafft.aln.fas
#cladeF.mafft.aln.fas
#cladeG.mafft.aln.fas
#cladeH.mafft.aln.fas
#cladeI.mafft.aln.fas
for(i in 1:9){
    fileName=paste("clade",cladeArray[i],".mafft.aln.fas",sep="")
    seq <- read.fas(file=fileName)
    #sic
    hsic<-haplotype(seq,indels="sic")
    dist<-hsic@d
    write.csv(dist,paste(fileName,".matrix.csv",sep=""))
    M <- mst(dist)
    write.csv(M,paste(fileName,".mst.csv",sep=""))
}
