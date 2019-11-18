#!/usr/bin/Rscript
##Rscript ~/script/hictools/1_select_loop.R  infile outpdf outsigloop pvalues_cut
##Rscript ~/script/hictools/1_select_loop.R loop.ihep2d_1_2.txt loop.ihep2d_1_2.pdf loop.ihep2d_1_2.sig_p0.01.txt 0.01


argv<-commandArgs(TRUE)
infile<-argv[1] ## must has header.Format:chr1,x1,x2,chr2,y1,y2,pets_num,ratio,back_pets
outpdf<-argv[2]
outsigloop<-argv[3]
pvalues_cut<-as.numeric(argv[4])


library(ggplot2)

SelectLoop<-function(infile, outpdf, outsigloop, pvalues_cut){
	data = read.table(infile, header = FALSE)
    colnames(data)=c("chr1","x1","x2","chr2","y1","y2","pets_num","ratio","back_pets")
	anchor_L_mid = data[,2] + (data[,3] - data[,2])/2
	anchor_R_mid = data[,5] + (data[,6] - data[,5])/2

	dis = anchor_R_mid - anchor_L_mid

	d1 = cbind(data, dis)

	m = mean(d1$dis)
	a=d1$pets_num*d1$dis/m
	max_pets_num=quantile(a,0.9999)
	a[which(a>max_pets_num)]=max_pets_num
	lam=mean(a)
    pvalues=ppois(a,lambda =lam, lower.tail=FALSE)
	pvalues=pvalues*length(pvalues)
	p_log = -log10(pvalues)
	v1 = cbind(data, pvalues, p_log)
	v1$pvalues[which(v1$pvalues>50)]=50

	sig=v1[which(v1$pvalues<=pvalues_cut),]
	print(length(sig[,1]))
	write.table(sig,file=outsigloop,sep="\t",quote=FALSE,row.names = FALSE)

	p = ggplot(v1,aes(x=dis,y=pets_num,colour=p_log)) + geom_point() + ylim(0,500) + xlim(0,1000000) + scale_color_gradient(low="green",high="red") + theme_classic()
	p = p + labs(title= "loop Pvalue of pets_num/distance ", x = "distance", y = "-log10(P)")
	p = p + theme_classic()
	p = p + theme(axis.text.x=element_text(size=rel(1.5), color = "black"),
	        	  axis.text.y=element_text(size=rel(1.5), color = "black"),
	        	  axis.title.x=element_text(size=rel(1.5), color = "black"),
	        	  axis.title.y=element_text(size=rel(1.5), color = "black"),
	        	  plot.title=element_text(size=rel(2),hjust=0.1, color = "black"),
	        	  panel.border=element_rect(colour="black",fill=NA,size=1))
	pdf(outpdf,width=6.5,height=5)
	print(p)
	dev.off()
}


SelectLoop(infile, outpdf, outsigloop, pvalues_cut)

