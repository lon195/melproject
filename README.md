##Assignment 3:  Start by identifying each section in the code below.  Then annotate as many lines as you can using the 
##comment symbol (#).  If you need help, try googling your question.  Also use the following references:  

###https://www.statmethods.net/r-tutorial/index.html
###https://rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf
###https://rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf
##Don't worry about sections 5 and 6 for now.
##Also - to make the P value comparison you will need to install the package ggpubr.  Let me know if you have issues.


##loading programs####
library(cowplot) # add on to ggplot
library(ggplot2) # creating a new ggplot
library(tidyverse) # allows data subsetting and transformation
install.packages("ggpubr") # allows customization of ggplot
library(ggpubr) 
install.packages("devtools")
devtools::install_github("username/packagename")

#importdataset-->session-->workingdirectory-->choose directory
read.csv(data2) # asks R to run function

##defining functions####
se <- function(x) sqrt(var(x, na.rm=TRUE)/length(x)) # defining std error and skipping over any NA values with na.rm=TRUE


data_summary<-function(x){
  m<-median(x)
  ymin<-m-se(x)
  ymax<-m+se(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
} # function for summary stats for ploting one stdev

data_summary1<-function(x){
  m<-median(x)
  ymin<-m-(IQR(x)/2)
  ymax<-m+(IQR(x)/2)
  return(c(y=m, ymin=ymin,ymax=ymax))
} # function for summary stats for ploting range of data

data_summary2<-function(x){
  m<-mean(x)
  ymin<-m-se(x)
  ymax<-m+se(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}# function for summary stats for ploting range of standard error


# defining data frames 
analyzer<-function(data,gene,control,date,plot_arg,outdir){
  if (dir.exists(outdir)==F) {dir.create(outdir)} # checks to see if the output directory exits.  If not, it creates the output directory
  df<-read.csv(data, header = TRUE)# reads the data file identified in the function and assigns it to the dataframe, df
  df_1<-df[which(df$censor=="n"),]  # subsets the dataframe to include non-censored observations and assigns to df_1
  df_2<-df_1[which(df_1$gene_tested==gene),] # subsets the dataframe df_1 to include observations where the variable gene_tested matches the function parameter gene.
  df_3<-df_2[which(df_2$experiment_date==date),] # subsets the dataframe df_2 to include observations where the variable gene-tested occurs on the date specified and assigns it as df_3
  df_4c<-df_3[which(df_3$condition==control),] # subsets the dataframe df_3 to include observations where the variable gene_tested doesn't match the gene tested as the control and assigns it as df_4c
  df_4e<-df_3[which(df_3$condition==gene),] # subsets the dataframe df_3 to include observations where the variable gene_tested matches the gene tested as the gene and assigns it as df_4e
  sink(paste0(outdir,"/stat_report.txt")) # sink to store the output in a file
  cat(paste0("Stat Report for ",gene," Vs ",control," on ",date,"\n\n")) #cat combines multiple items into a continuous output
  print(df_3 %>% group_by(condition) %>% summarise(n=n(), mean = mean(blinded_score), sd = sd(blinded_score), se = se(blinded_score), median = median(blinded_score), IQR = IQR(blinded_score))) # group blinded score by condition
  con_bs<-df_4c$blinded_score #subsets the dataframe to include observations (blinded scores) of the control assigned as con_bs
  exp_bs<-df_4e$blinded_score #subsets the dataframe to include observations (blinded scores) of the gene assigned as exp_bs
  if (length(con_bs)>1 && length(exp_bs>1)) {# condition blinded scores must be >1; Kolmogorov-Smirnov Tests; not paired b/c means not normally distributed
    print(ks.test(con_bs,exp_bs))
    print(t_equal<-t.test(con_bs,exp_bs,alternative="two.sided",paired=FALSE,var.equal=TRUE)) # var.equal = TRUE will estimate variance using pooled variance
    print(t_unequal<-t.test(con_bs,exp_bs,alternative="two.sided",paired=FALSE,var.equal=FALSE)) # var.equal = FALSE assumes inequal variance between datasets
    print(wilcox<-wilcox.test(con_bs,exp_bs,alternative="two.sided"))
  }#null hypothesis: distribution function of con_bs = distribution function of exp_bs
  sink() # diverts R output and stops it
  df_3$condition<-relevel(df_3$condition,control) # converts df for plotting
  plot<-ggplot(data=df_3, aes(x=condition, y=blinded_score, fill = condition)) # identifies which variable should be mapped to which aesthetics (aes)
  plot<-plot+
    geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03)) + # stacking direction of dots=stackdir along Y-axis
    scale_fill_manual(values=c("#3C5488","#DC0000")) + 
    scale_color_manual(values=c("black","black"))+
    stat_summary(fun.data=data_summary2, color="black", size=0.5, width=0.3, alpha=0.3,geom="crossbar")+ # error bars
    stat_compare_means(label = "p.signif", comparisons = list(c(control,gene)))+ # add p-values
    theme_cowplot(font_size = 11)+
    labs(y="Blinded Score", x=NULL)+ # labels no X Y=blinded score
    ylim(c(0,1.1*max(df_3$blinded_score)))+ #y limit based on values for df_3 blinded scores
    theme(legend.position = "none") 
  if (plot_arg==TRUE){
    save_plot(plot = plot, filename=paste0(outdir,"/",gene,".pdf"),base_width = 3, base_height = 2.7)
  }
  #return(df_3)
}

# plotting data for dlx4a vs control on 1 day
analyzer(data = "data2.csv", gene = "dlx4a", control = "mcs",date = 20181102, plot_arg = TRUE,outdir = "dlx4a")

# delay re-encoding of strings
toc<-read.csv("data2.csv", header = TRUE, stringsAsFactors = F) # import csv and string as factors imports column names as character data (false) or as factors (true)
toc<-toc %>% select(experiment_date,gene_tested) %>% distinct() # select certain rows i.e. same data and gene tested

# plotting data for each gene on each day
for (i in 1:nrow(toc)) {
  analyzer(data = "data2.csv", gene = toc[i,2], control = "mcs", date = toc[i,1], plot_arg = T, outdir = paste0(toc[i,2],"_",toc[i,1]))
} # specify with nrow to only satisy the conditions of experiment date and gene tested

