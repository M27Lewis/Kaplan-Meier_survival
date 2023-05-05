library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(GSVA)
library(survival)
library(rms)

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

# Assign experiment name
experiment <- "GSE96058_KM_analysis"

# Create directories to save output files
plots_path <- file.path(paste0("plots/", experiment))

if(!dir.exists(plots_path)){
  dir.create(plots_path)
}

tables_path <- file.path(paste0("tables/", experiment))

if(!dir.exists(tables_path)){
  dir.create(tables_path)
}

# Import gene expression and clinical data for KM analysis
sb_genez <- readRDS("data/SCAN-B_GSE96058/SB_subset_gene_expression_z-scored.rds")

sb_clin <- readRDS("data/SCAN-B_GSE96058/SB_reduced_clinical_data.rds")

sb_clin <- sb_clin |> 
  rename(PATIENT_ID = Sample_title)

#Import gene lists and prepare for analysis
gene_path <- "data/gene_signatures/SE_DEGs_clean.txt" #complete list of DEGs from CRISPRi Screen 2.0 experiment

genes <- read.delim(gene_path, header = TRUE, sep = "\t") |> 
  dplyr::rename("cohort" = "SE") |> 
  dplyr::rename("gene" = "X") |> 
  filter(!cohort %in% c("CSF1", "CSF1e")) |>  #filter to remove positive controls
  filter(baseMean > 10) |>  #remove genes with basemean expression < 10 (same as RNA-seq analysis)
  filter(log2FoldChange < 0) |> #remove genes that increase in expression when SE is knocked down
  filter(padj < 0.05) #only DEG adj pval < 0.05

#Iterate through gene list to assign median score for each patient
cohort_list <- unique(genes$cohort)

results <- as.data.frame(matrix(
  ncol = length(unique(genes$cohort)), 
  nrow = nrow(sb_clin))) #build empty df for results, # of rows matches # of patients, # of columns matches # of cohorts

colnames(results) <- cohort_list #Assign SE names to the column names of results DF

results$PATIENT_ID <- colnames(sb_genez) #Create CLID column to keep track of patient IDs



for (i in 1:length(cohort_list)){
  gene_list <- genes |> 
    filter(cohort == cohort_list[i]) #subset to rows that contain genes listed in the [i] group
  surv_data <- as.data.frame(t(sb_genez |> 
                                 filter(rownames(sb_genez) %in% gene_list$gene))) 
  
  results[,i] <- apply(surv_data, 1, median)
  
}

results[,1:(ncol(results)-1)] <- scale(results[,1:(ncol(results)-1)]) #Standardize the median results by column (mean=0, sd=1), but skip PATIENT_ID column

x <- merge(sb_clin, results, by = "PATIENT_ID")
