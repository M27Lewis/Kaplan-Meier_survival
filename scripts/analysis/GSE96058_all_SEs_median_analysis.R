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


# Convert survival days to years for KM analysis
x <- x |> 
  mutate(DRFi_days = DRFi_days/365.2425) |> 
  rename(DRFi_years = DRFi_days) |> 
  mutate(OS_days = OS_days/365.2425) |> 
  rename(OS_years = OS_days) |> 
  mutate(RFi_days = RFi_days/365.2425) |> 
  rename(RFi_years = RFi_days) |> 
  mutate(BCFi_days = BCFi_days/365.2425) |> 
  rename(BCFi_years = BCFi_days)


#Run KM analysis and save plots using DRFi
this.name="All_SEs_padj_05"
pdf(paste(plots_path, "/", this.name, "_Median_scaled_in_SCAN-B_DRFi_KM.pdf",sep=""), width=12, height=7) #Default is W = 12 and H = 7 for larger PDFs
par(mgp=c(1.3,.35,.0),mai=c(.5,.5,.4,.2),mfrow=c(1,2))



for(i in 21:dim(x)[2])
{
  x<-subset(x,x[,14]!="NA") #removing NAs from the survival column
  dim(x)
  x<-subset(x,x[,13]!="#VALUE!") #removing non-numeric values from survival time
  dim(x)
  #x<-subset(x,x[,i]!="NA") #removing any general NAs
  #dim(x)
  
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
         x$DRFi_event, 
         x$DRFi_years, 
         ylab="DRFi (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$DRFi_years, x$DRFi_event))[[1]],3)	
  fit<-summary(coxph(Surv(x$DRFi_years, x$DRFi_event)~x[,i]))	
  if(fit$coef[,5]<0.05)
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n", text.col="red")	
  }else
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n")	
  }
  
  kmPlot(G.3, 
         x$DRFi_event, 
         x$DRFi_years, 
         ylab="DRFi (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","blue","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$DRFi_years, x$DRFi_event))[[1]],3)	
  fit<-summary(coxph(Surv(x$DRFi_years, x$DRFi_event)~x[,i]))	
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

for(i in 21:dim(x)[2])
{
  x<-subset(x,x[,14]!="NA") #removing NAs from the survival column
  dim(x)
  x<-subset(x,x[,13]!="#VALUE!") #removing non-numeric values from survival time
  dim(x)
  #x<-subset(x,x[,i]!="NA") #removing any general NAs
  #dim(x)
  
  fit<-summary(coxph(Surv(x$DRFi_years, x$DRFi_event)~x[,i]))	
  
  HR_df[i-20, "SE"] <- colnames(x[i])
  HR_df[i-20,"Cox_hazard_ratio"] <- signif(fit$coef[,2],4)
  HR_df[i-20,"95CI_min"] <- signif(fit$conf.int[,3],3)
  HR_df[i-20,"95CI_max"] <- signif(fit$conf.int[,4],3)
  HR_df[i-20,"Cox_p"] <- signif(fit$coef[,5],3)
  
}

write_tsv(as.data.frame(HR_df), paste(tables_path, "/Cox_hazard_ratio_results_all_SEs_padj_0.05_Median_scaled_in_SCAN-B_DRFi_KM.txt", sep = ""))


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


ggsave(paste(plots_path, "/All_SEs_padj_0.05_Forest_plot_SCAN-B_DRFi.pdf", sep = ""), plot = fp1, width = 6, height = 5.5, units = "in", dpi = 600)
ggsave(paste(plots_path, "/All_SEs_padj_0.05_Forest_plot_SCAN-B_DRFi.png", sep = ""), plot = fp1, width = 6, height = 5.5, units = "in", dpi = 600)



## Run KM analysis and save plots using OS
this.name="All_SEs_padj_05"
pdf(paste(plots_path, "/", this.name, "_Median_scaled_in_SCAN-B_OS_KM.pdf",sep=""), width=12, height=7) #Default is W = 12 and H = 7 for larger PDFs
par(mgp=c(1.3,.35,.0),mai=c(.5,.5,.4,.2),mfrow=c(1,2))



for(i in 21:dim(x)[2])
{
  x<-subset(x,x[,16]!="NA") #removing NAs from the survival column
  dim(x)
  x<-subset(x,x[,15]!="#VALUE!") #removing non-numeric values from survival time
  dim(x)
  #x<-subset(x,x[,i]!="NA") #removing any general NAs
  #dim(x)
  
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
         x$OS_event, 
         x$OS_years, 
         ylab="OS (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$OS_years, x$OS_event))[[1]],3)	
  fit<-summary(coxph(Surv(x$OS_years, x$OS_event)~x[,i]))	
  if(fit$coef[,5]<0.05)
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n", text.col="red")	
  }else
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n")	
  }
  
  kmPlot(G.3, 
         x$OS_event, 
         x$OS_years, 
         ylab="OS (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","blue","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$OS_years, x$OS_event))[[1]],3)	
  fit<-summary(coxph(Surv(x$OS_years, x$OS_event)~x[,i]))	
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
HR_df2 <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("SE", "Cox_hazard_ratio", "95CI_min", "95CI_max", "Cox_p"))

for(i in 21:dim(x)[2])
{
  x<-subset(x,x[,16]!="NA") #removing NAs from the survival column
  dim(x)
  x<-subset(x,x[,15]!="#VALUE!") #removing non-numeric values from survival time
  dim(x)
  #x<-subset(x,x[,i]!="NA") #removing any general NAs
  #dim(x)
  
  fit<-summary(coxph(Surv(x$OS_years, x$OS_event)~x[,i]))	
  
  HR_df2[i-20, "SE"] <- colnames(x[i])
  HR_df2[i-20,"Cox_hazard_ratio"] <- signif(fit$coef[,2],4)
  HR_df2[i-20,"95CI_min"] <- signif(fit$conf.int[,3],3)
  HR_df2[i-20,"95CI_max"] <- signif(fit$conf.int[,4],3)
  HR_df2[i-20,"Cox_p"] <- signif(fit$coef[,5],3)
  
}

write_tsv(as.data.frame(HR_df2), paste(tables_path, "/Cox_hazard_ratio_results_all_SEs_padj_0.05_Median_scaled_in_SCAN-B_OS_KM.txt", sep = ""))


# Forest plot of hazard ratios for all SEs
HR_df2$SE <- factor(HR_df2$SE, levels = HR_df2$SE[order(HR_df2$Cox_hazard_ratio)]) #Reorder so that largest HR is on top and goes in descending order

fp2 <- ggplot(HR_df2, mapping = aes(x = Cox_hazard_ratio, y = SE)) +
  geom_point(colour = "black", shape = 10, size = 4) + #Selecting point color, shape, and size
  geom_errorbarh(aes(xmin = `95CI_min`, xmax = `95CI_max`), height = 0.5) + #Creates the lower and upper bars
  geom_vline(xintercept = 1, linetype = "dashed") + #Adds intercept line at HR = 1
  geom_text(aes(x = (min(`95CI_min`) - 0.2), label = Cox_p), size = 3) + #Prints the cox pval on the plot to the left of the smallest value
  scale_x_continuous(expand = expansion(mult = .2)) + #Auto scales the plot size with a 10% buffer
  labs(x = 'Hazard Ratio', y = '') + #Setting the X and Y labels
  theme(legend.position = 'none') + #Removing the legend
  theme_Publication()

fp2


ggsave(paste(plots_path, "/All_SEs_padj_0.05_Forest_plot_SCAN-B_OS.pdf", sep = ""), plot = fp2, width = 6, height = 5.5, units = "in", dpi = 600)
ggsave(paste(plots_path, "/All_SEs_padj_0.05_Forest_plot_SCAN-B_OS.png", sep = ""), plot = fp2, width = 6, height = 5.5, units = "in", dpi = 600)


## Run KM analysis and save plots using RFi
this.name="All_SEs_padj_05"
pdf(paste(plots_path, "/", this.name, "_Median_scaled_in_SCAN-B_RFi_KM.pdf",sep=""), width=12, height=7) #Default is W = 12 and H = 7 for larger PDFs
par(mgp=c(1.3,.35,.0),mai=c(.5,.5,.4,.2),mfrow=c(1,2))



for(i in 21:dim(x)[2])
{
  x<-subset(x,x[,18]!="NA") #removing NAs from the survival column
  dim(x)
  x<-subset(x,x[,17]!="#VALUE!") #removing non-numeric values from survival time
  dim(x)
  #x<-subset(x,x[,i]!="NA") #removing any general NAs
  #dim(x)
  
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
         x$RFi_event, 
         x$RFi_years, 
         ylab="RFi (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$RFi_years, x$RFi_event))[[1]],3)	
  fit<-summary(coxph(Surv(x$RFi_years, x$RFi_event)~x[,i]))	
  if(fit$coef[,5]<0.05)
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n", text.col="red")	
  }else
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n")	
  }
  
  kmPlot(G.3, 
         x$RFi_event, 
         x$RFi_years, 
         ylab="RFi (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","blue","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$RFi_years, x$RFi_event))[[1]],3)	
  fit<-summary(coxph(Surv(x$RFi_years, x$RFi_event)~x[,i]))	
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
HR_df3 <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("SE", "Cox_hazard_ratio", "95CI_min", "95CI_max", "Cox_p"))

for(i in 21:dim(x)[2])
{
  x<-subset(x,x[,18]!="NA") #removing NAs from the survival column
  dim(x)
  x<-subset(x,x[,17]!="#VALUE!") #removing non-numeric values from survival time
  dim(x)
  #x<-subset(x,x[,i]!="NA") #removing any general NAs
  #dim(x)
  
  fit<-summary(coxph(Surv(x$RFi_years, x$RFi_event)~x[,i]))	
  
  HR_df3[i-20, "SE"] <- colnames(x[i])
  HR_df3[i-20,"Cox_hazard_ratio"] <- signif(fit$coef[,2],4)
  HR_df3[i-20,"95CI_min"] <- signif(fit$conf.int[,3],3)
  HR_df3[i-20,"95CI_max"] <- signif(fit$conf.int[,4],3)
  HR_df3[i-20,"Cox_p"] <- signif(fit$coef[,5],3)
  
}

write_tsv(as.data.frame(HR_df3), paste(tables_path, "/Cox_hazard_ratio_results_all_SEs_padj_0.05_Median_scaled_in_SCAN-B_RFi_KM.txt", sep = ""))


# Forest plot of hazard ratios for all SEs
HR_df3$SE <- factor(HR_df3$SE, levels = HR_df3$SE[order(HR_df3$Cox_hazard_ratio)]) #Reorder so that largest HR is on top and goes in descending order

fp3 <- ggplot(HR_df3, mapping = aes(x = Cox_hazard_ratio, y = SE)) +
  geom_point(colour = "black", shape = 10, size = 4) + #Selecting point color, shape, and size
  geom_errorbarh(aes(xmin = `95CI_min`, xmax = `95CI_max`), height = 0.5) + #Creates the lower and upper bars
  geom_vline(xintercept = 1, linetype = "dashed") + #Adds intercept line at HR = 1
  geom_text(aes(x = (min(`95CI_min`) - 0.2), label = Cox_p), size = 3) + #Prints the cox pval on the plot to the left of the smallest value
  scale_x_continuous(expand = expansion(mult = .2)) + #Auto scales the plot size with a 10% buffer
  labs(x = 'Hazard Ratio', y = '') + #Setting the X and Y labels
  theme(legend.position = 'none') + #Removing the legend
  theme_Publication()

fp3


ggsave(paste(plots_path, "/All_SEs_padj_0.05_Forest_plot_SCAN-B_RFi.pdf", sep = ""), plot = fp3, width = 6, height = 5.5, units = "in", dpi = 600)
ggsave(paste(plots_path, "/All_SEs_padj_0.05_Forest_plot_SCAN-B_RFi.png", sep = ""), plot = fp3, width = 6, height = 5.5, units = "in", dpi = 600)


## Run KM analysis and save plots using BCFi
this.name="All_SEs_padj_05"
pdf(paste(plots_path, "/", this.name, "_Median_scaled_in_SCAN-B_BCFi_KM.pdf",sep=""), width=12, height=7) #Default is W = 12 and H = 7 for larger PDFs
par(mgp=c(1.3,.35,.0),mai=c(.5,.5,.4,.2),mfrow=c(1,2))



for(i in 21:dim(x)[2])
{
  x<-subset(x,x[,20]!="NA") #removing NAs from the survival column
  dim(x)
  x<-subset(x,x[,19]!="#VALUE!") #removing non-numeric values from survival time
  dim(x)
  #x<-subset(x,x[,i]!="NA") #removing any general NAs
  #dim(x)
  
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
         x$BCFi_event, 
         x$BCFi_years, 
         ylab="BCFi (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$BCFi_years, x$BCFi_event))[[1]],3)	
  fit<-summary(coxph(Surv(x$BCFi_years, x$BCFi_event)~x[,i]))	
  if(fit$coef[,5]<0.05)
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n", text.col="red")	
  }else
  {
    legend("topright",legend=paste("Cox HR = ",signif(fit$coef[,2],4), "\n95% CI (", signif(fit$conf.int[,3],3), ", ", signif(fit$conf.int[,4],3), ")", "\nCox P = ", signif(fit$coef[,5],3),"\nC-index = ",cindex), bty="n")	
  }
  
  kmPlot(G.3, 
         x$BCFi_event, 
         x$BCFi_years, 
         ylab="BCFi (Probability)",
         xlab="Years",
         lwid=1,
         lineColors=c("green","blue","red"), 
         varName=names(x)[i])
  
  cindex<- signif(1-rcorr.cens(x[,i],Surv(x$BCFi_years, x$BCFi_event))[[1]],3)	
  fit<-summary(coxph(Surv(x$BCFi_years, x$BCFi_event)~x[,i]))	
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
HR_df4 <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("SE", "Cox_hazard_ratio", "95CI_min", "95CI_max", "Cox_p"))

for(i in 21:dim(x)[2])
{
  x<-subset(x,x[,20]!="NA") #removing NAs from the survival column
  dim(x)
  x<-subset(x,x[,19]!="#VALUE!") #removing non-numeric values from survival time
  dim(x)
  #x<-subset(x,x[,i]!="NA") #removing any general NAs
  #dim(x)
  
  fit<-summary(coxph(Surv(x$BCFi_years, x$BCFi_event)~x[,i]))	
  
  HR_df4[i-20, "SE"] <- colnames(x[i])
  HR_df4[i-20,"Cox_hazard_ratio"] <- signif(fit$coef[,2],4)
  HR_df4[i-20,"95CI_min"] <- signif(fit$conf.int[,3],3)
  HR_df4[i-20,"95CI_max"] <- signif(fit$conf.int[,4],3)
  HR_df4[i-20,"Cox_p"] <- signif(fit$coef[,5],3)
  
}

write_tsv(as.data.frame(HR_df4), paste(tables_path, "/Cox_hazard_ratio_results_all_SEs_padj_0.05_Median_scaled_in_SCAN-B_BCFi_KM.txt", sep = ""))


# Forest plot of hazard ratios for all SEs
HR_df4$SE <- factor(HR_df4$SE, levels = HR_df4$SE[order(HR_df4$Cox_hazard_ratio)]) #Reorder so that largest HR is on top and goes in descending order

fp4 <- ggplot(HR_df4, mapping = aes(x = Cox_hazard_ratio, y = SE)) +
  geom_point(colour = "black", shape = 10, size = 4) + #Selecting point color, shape, and size
  geom_errorbarh(aes(xmin = `95CI_min`, xmax = `95CI_max`), height = 0.5) + #Creates the lower and upper bars
  geom_vline(xintercept = 1, linetype = "dashed") + #Adds intercept line at HR = 1
  geom_text(aes(x = (min(`95CI_min`) - 0.2), label = Cox_p), size = 3) + #Prints the cox pval on the plot to the left of the smallest value
  scale_x_continuous(expand = expansion(mult = .2)) + #Auto scales the plot size with a 10% buffer
  labs(x = 'Hazard Ratio', y = '') + #Setting the X and Y labels
  theme(legend.position = 'none') + #Removing the legend
  theme_Publication()

fp4


ggsave(paste(plots_path, "/All_SEs_padj_0.05_Forest_plot_SCAN-B_BCFi.pdf", sep = ""), plot = fp4, width = 6, height = 5.5, units = "in", dpi = 600)
ggsave(paste(plots_path, "/All_SEs_padj_0.05_Forest_plot_SCAN-B_BCFi.png", sep = ""), plot = fp4, width = 6, height = 5.5, units = "in", dpi = 600)