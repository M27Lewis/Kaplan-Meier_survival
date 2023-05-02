library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(GSVA)
library(survival)
library(rms)
library('org.Hs.eg.db')

Surv<- survival::Surv

###########################################################
##########  Function Definition               #############
###########################################################

kmPlot<-function(x,
                 event,
                 stime,
                 varName="",
                 ymin=0,
                 lineColors=NA,
                 nclasses=NA,
                 ci=F,
                 ylab="OS (Probability)",
                 xlab="Years",
                 plegloc="bottomright",
                 lwid=1,
                 citime=c(5,10)){
  
  mainLabel <- varName
  event<-as.numeric(as.vector(event))
  stime<-as.numeric(as.vector(stime))
  if(is.numeric(x)){
    tclass <- factor(x[!is.na(x)])
    event <- event[!is.na(x)]
    stime <- stime[!is.na(x)]
  }else{
    tclass <- factor(x[x!=""])
    event <- event[x!=""]
    stime <- stime[x!=""]
  }
  
  if(is.na(nclasses)){
    nclasses<-nlevels(tclass)
  }
  
  if(length(lineColors)<=1){
    lineColors<-seq(1,nclasses)
  }
  
  y<-survfit(Surv(stime, event)~tclass) #applying the survival function
  
  plot(y,cex.main=0.7,
       col=lineColors,
       main=mainLabel,
       ylim=c(ymin,1),
       ylab=ylab,
       xlab=xlab,
       lwd=lwid)
  
  sy<-summary(y)
  
  if(ci==T)
  {
    for(j in 2:(nclasses+1))
    {
      thislevel<-which(sy$strata==levels(sy$strata)[j-1])
      if(levels(sy$strata)[j-1] != "Basal")
      {
        for(k in 1:length(citime))
        {
          closest<-min(abs(citime[k]-sy$time[thislevel]))==abs(citime[k]-sy$time[thislevel])
          low<-sy$lower[thislevel][closest]
          high<-sy$upper[thislevel][closest]
          segments(citime[k],low,citime[k],high,lty=2)
          segments(citime[k]-0.3,low,citime[k]+0.3,low)
          segments(citime[k]-0.3,high,citime[k]+0.3,high)
          print(paste(levels(sy$strata)[j-1],k,low,sy$surv[thislevel][closest],high))
        }
      }
    }
  }
  legendtxt<- paste(levels(tclass)," ",table(tclass,event)[,2],"/",y$n,sep="")
  legend(0,0.3,legend=legendtxt,col=lineColors,bty="n",lty=rep(1,nclasses),lwd=lwid)
  
  pvalue<-1-pchisq(survdiff(Surv(stime, event)~tclass)$chisq,nclasses-1)
  if(pvalue<0.05)
  {
    legend(plegloc,legend=paste("log rank p = ",signif(pvalue,3),sep=""),cex=1.2,text.col="red",bty="n")
  }else{
    legend(plegloc,legend=paste("log rank p=",signif(pvalue,3),sep=""),cex=1.2,bty="n")
  }	
}

#####################################

# Custom ggplot theme
theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, margin = margin(b = 5)),
            text = element_text(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = margin(0, 0, 0, 0),
            legend.title = element_text(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}
##############################


# Import datasets for analysis
MD <- readRDS("data/GSE25066/MDACC_Clinical_Expression_data.rds")

MDF <- MD |> 
  dplyr::slice(46:n()) |> #Remove first 44 rows
  filter(!grepl("///", CLID)) |> #Remove lines with multiple probe hits
  mutate_all(function(x) as.numeric(as.character(x))) #Convert all to numeric

rownames(MDF) <- MDF$CLID

MDF <- MDF[,-1] #Remove CLID Row containing Entrez IDs

#Importing the reduced survival data from Chris to append my gsva SE signature scores to
MDA_data <- read.delim("data/GSE25066/survival_data_CF.txt", sep="\t",header=T)

#Remove punctuation to make CLIDs match
MDA_data$CLID <- gsub('[[:punct:]]','', MDA_data$CLID)


#Import gene lists and prepare for analysis
gene_path <- "data/gene_signatures/SE_DEGs_clean.txt" #complete list of DEGs from CRISPRi Screen 2.0 experiment

genes <- read.delim(gene_path, header = TRUE, sep = "\t") |> 
  dplyr::rename("cohort" = "SE") |> 
  dplyr::rename("gene" = "X") |> 
  filter(!cohort %in% c("CSF1", "CSF1e")) |>  #filter to remove positive controls
  filter(baseMean > 10) |>  #remove genes with basemean expression < 10 (same as RNA-seq analysis)
  filter(log2FoldChange < 0) |> #remove genes that increase in expression when SE is knocked down
  filter(padj < 0.05) #only DEG adj pval < 0.05


#For loop to iterate through signature list and find medians
cohort_list <- unique(genes$cohort)

results <- as.data.frame(matrix(
  ncol = length(unique(genes$cohort)), 
  nrow = ncol(MDF))) #build empty df for results, # of rows matches # of patients, # of columns matches # of cohorts

colnames(results) <- cohort_list #Assign SE names to the column names of results DF

results$CLID <- colnames(MDF) #Create CLID column to keep track of patient IDs

genes$Entrez <- mapIds(org.Hs.eg.db, genes$gene, 'ENTREZID', 'SYMBOL') #Create a column with the Entrez ID that matches each gene (for appending survival data later)


for (i in 1:length(cohort_list)){
  gene_list <- genes |> 
    filter(cohort == cohort_list[i])
  surv_data <- as.data.frame(t(MDF |> 
                                 filter(rownames(MDF) %in% gene_list$Entrez)))
  
  results[,i] <- apply(surv_data, 1, median)
  
}

results[,1:(ncol(results)-1)] <- scale(results[,1:(ncol(results)-1)]) #Standardize the median results by column (mean=0, sd=1), but skip CLID column

results$CLID <- gsub('[[:punct:]]','', results$CLID)

x <- merge(MDA_data, results, by = "CLID")


# Create and save KM plots
this.name="All_SEs_GSE25066_padj_0.05"
par(mgp=c(1.3,.35,.0),mai=c(.5,.5,.4,.2))
pdf(paste("plots/", this.name, "_Median_scaled_in_MDACC_DRFS_KM.pdf",sep=""), width=12, height=7)

par(mfrow=c(1,2))
for (i in 8:dim(x)[2])
{
  
  quantile.00<-as.numeric(quantile(x[,i], 0.00))
  quantile.25<-as.numeric(quantile(x[,i], 0.25))
  quantile.33<-as.numeric(quantile(x[,i], 0.3333))
  quantile.50<-as.numeric(quantile(x[,i], 0.50))
  quantile.66<-as.numeric(quantile(x[,i], 0.6666))
  quantile.75<-as.numeric(quantile(x[,i], 0.75))
  quantile.100<-as.numeric(quantile(x[,i], 0.100))
  
  
  G.2<- cut(x[,i],c(-1000,quantile.50,1000),c("low","high"))
  G.2<- cut(x[,i],c(-1000,quantile.50,1000),c("low","high"))
  G.3<- cut(x[,i],c(-1000,quantile.33,quantile.66,1000),c("low","med","high"))
  G.4<- cut(x[,i],c(-1000,quantile.25,quantile.50,quantile.75,1000),c("low","med.low","med.high","high"))
  
  kmPlot(G.2, 
         x$drfs1event0censored, 
         x$drfseventimeyears, 
         ylab="DRFS (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$drfseventimeyears, x$drfs1event0censored))[[1]],3)	
  fit<-summary(coxph(Surv(x$drfseventimeyears, x$drfs1event0censored)~x[,i]))	
  if(fit$coef[,5]<0.05)
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n", text.col="red")	
  }else
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n")	
  }
  
  kmPlot(G.3, 
         x$drfs1event0censored, 
         x$drfseventimeyears, 
         ylab="DRFS (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","blue","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$drfseventimeyears, x$drfs1event0censored))[[1]],3)	
  fit<-summary(coxph(Surv(x$drfseventimeyears, x$drfs1event0censored)~x[,i]))	
  #legend("topright",legend=paste("C-index = ",cindex,  "\nCox HR = ",signif(fit$coef[,2],4), "\nCox p = ", signif(fit$coef[,5],3)), bty="n")	
  
  if(fit$coef[,5]<0.05)
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n", text.col="red")	
  }else
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n")	
  }
  
}


dev.off()