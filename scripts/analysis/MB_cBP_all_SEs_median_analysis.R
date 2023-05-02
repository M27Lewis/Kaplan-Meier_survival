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


#Import datasets for analysis
mb_cbp_z <- readRDS("data/METABRIC_cBioPortal/METABRIC_cBioPortal_z-scored_data.rds")

mb_cbp_clin <- readRDS("data/METABRIC_cBioPortal/METABRIC_cBioPortal_clinical_data_processed.rds")

#Modify datasets to ensure matching patient IDs
ids <- unique(mb_cbp_clin$PATIENT_ID) #clinical patient IDs

ids2 <- unique(colnames(mb_cbp_z)) #expression patient IDs

c(setdiff(ids2, ids), setdiff(ids, ids2)) #patient ID(s) that do not appear in both datasets

mb_cbp_z <- mb_cbp_z |> 
  select(!"MB-5130") #Removing this ID since it does not appear in clinical data (might have NAs)

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
  nrow = nrow(mb_cbp_clin))) #build empty df for results, # of rows matches # of patients, # of columns matches # of cohorts

colnames(results) <- cohort_list #Assign SE names to the column names of results DF

results$PATIENT_ID <- colnames(mb_cbp_z) #Create CLID column to keep track of patient IDs



for (i in 1:length(cohort_list)){
  gene_list <- genes |> 
    filter(cohort == cohort_list[i]) #subset to rows that contain genes listed in the [i] group
  surv_data <- as.data.frame(t(mb_cbp_z |> 
                                 filter(rownames(mb_cbp_z) %in% gene_list$gene))) 
  
  results[,i] <- apply(surv_data, 1, median)
  
}

results[,1:(ncol(results)-1)] <- scale(results[,1:(ncol(results)-1)]) #Standardize the median results by column (mean=0, sd=1), but skip PATIENT_ID column

x <- merge(mb_cbp_clin, results, by = "PATIENT_ID")


#Adjust survival data to correct class and cap at 10 years of survival data
x$OS_MONTHS <- as.numeric(x$OS_MONTHS) #convert column to numeric

x$OS_MONTHS <- x$OS_MONTHS / 12


x$OS_STATUS <- gsub("\\:.*", "", x$OS_STATUS) #remove the colon and everything after

x$OS_STATUS <- as.numeric(x$OS_STATUS) #convert to numeric

x$OS_MONTHS <- ifelse(x$OS_MONTHS >= 10, 10, x$OS_MONTHS) #cap survival time at 10 years

x <- x |> 
  mutate(OS_STATUS = ifelse(OS_MONTHS == 10 & OS_STATUS == 1, 0, OS_STATUS))



#Run KM analysis and save plots
this.name="All_SEs_padj_05"
pdf(paste("plots/", Sys.Date(), "_", this.name, "_Median_scaled_in_METABRIC_cBioPortal_OS_KM.pdf",sep=""), width=12, height=7) #Default is W = 12 and H = 7 for larger PDFs
par(mgp=c(1.3,.35,.0),mai=c(.5,.5,.4,.2),mfrow=c(1,2))



for(i in 8:dim(x)[2])
{
  x<-subset(x,x[,3]!="NA") #removing NAs from the survival column
  dim(x)
  x<-subset(x,x[,2]!="#VALUE!") #removing non-numeric values from survival time
  dim(x)
  x<-subset(x,x[,i]!="NA") #removing any general NAs
  dim(x)
  
  quantile.00<-as.numeric(quantile(x[,i], 0.00))
  quantile.25<-as.numeric(quantile(x[,i], 0.25))
  quantile.33<-as.numeric(quantile(x[,i], 0.3333))
  quantile.50<-as.numeric(quantile(x[,i], 0.50))
  quantile.66<-as.numeric(quantile(x[,i], 0.6666))
  quantile.75<-as.numeric(quantile(x[,i], 0.75))
  quantile.100<-as.numeric(quantile(x[,i], 0.100))
  
  #Cutting the data by quantile, with the 1000s being the outer boundaries that will never be reached
  G.2<- cut(x[,i],c(-1000,quantile.50,1000),c("low","high"))
  G.2<- cut(x[,i],c(-1000,quantile.50,1000),c("low","high"))
  G.3<- cut(x[,i],c(-1000,quantile.33,quantile.66,1000),c("low","med","high"))
  G.4<- cut(x[,i],c(-1000,quantile.25,quantile.50,quantile.75,1000),c("low","med.low","med.high","high"))
  
  kmPlot(G.2, 
         x$OS_STATUS, 
         x$OS_MONTHS, 
         ylab="OS (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$OS_MONTHS, x$OS_STATUS))[[1]],3)	
  fit<-summary(coxph(Surv(x$OS_MONTHS, x$OS_STATUS)~x[,i]))	
  if(fit$coef[,5]<0.05)
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n", text.col="red")	
  }else
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n")	
  }
  
  kmPlot(G.3, 
         x$OS_STATUS, 
         x$OS_MONTHS, 
         ylab="OS (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","blue","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$OS_MONTHS, x$OS_STATUS))[[1]],3)	
  fit<-summary(coxph(Surv(x$OS_MONTHS, x$OS_STATUS)~x[,i]))	
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


# Calculate hazard ratios for each SE and write to table
HR_df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("SE", "Cox_hazard_ratio", "95CI_min", "95CI_max", "Cox_p"))

for(i in 8:dim(x)[2])
{
  x<-subset(x,x[,3]!="NA") #removing NAs from the survival column
  dim(x)
  x<-subset(x,x[,2]!="#VALUE!") #removing non-numeric values from survival time
  dim(x)
  x<-subset(x,x[,i]!="NA") #removing any general NAs
  dim(x)
  
  
  fit<-summary(coxph(Surv(x$OS_MONTHS, x$OS_STATUS)~x[,i]))
  
  HR_df[i-7, "SE"] <- colnames(x[i])
  HR_df[i-7,"Cox_hazard_ratio"] <- signif(fit$coef[,2],4)
  HR_df[i-7,"95CI_min"] <- signif(fit$conf.int[,3],3)
  HR_df[i-7,"95CI_max"] <- signif(fit$conf.int[,4],3)
  HR_df[i-7,"Cox_p"] <- signif(fit$coef[,5],3)
  
}

write_tsv(as.data.frame(HR_df), "tables/Cox_hazard_ratio_results_all_SEs_METABRIC.txt")


# Forest plot of hazard ratios for all SEs
HR_df$SE <- factor(HR_df$SE, levels = HR_df$SE[order(HR_df$Cox_hazard_ratio)]) #Reorder so that largest HR is on top and goes in descending order

fp1 <- ggplot(HR_df, mapping = aes(x = Cox_hazard_ratio, y = SE)) +
  geom_point(colour = "black", shape = 10, size = 4) + #Selecting point color, shape, and size
  geom_errorbarh(aes(xmin = `95CI_min`, xmax = `95CI_max`), height = 0.5) + #Creates the lower and upper bars
  geom_vline(xintercept = 1, linetype = "dashed") + #Adds intercept line at HR = 1
  geom_text(aes(x = (min(`95CI_min`) - 0.2), label = Cox_p), size = 3) + #Prints the cox pval on the plot to the left of the smallest value
  scale_x_continuous(expand = expansion(mult = .2)) + #Auto scales the plot size with a 10% buffer
  labs(x = 'Hazard Ratio', y = '') + #Setting the X and Y labels
  theme(legend.position = 'none') + #Removing the legend
  theme_Publication()

fp1


ggsave("plots/All_SEs_padj_0.05_Forest_plot_METABRIC.pdf", plot = fp1, width = 6, height = 5.5, units = "in", dpi = 600)
ggsave("plots/All_SEs_padj_0.05_Forest_plot_METABRIC.png", plot = fp1, width = 6, height = 5.5, units = "in", dpi = 600)
