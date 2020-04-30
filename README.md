##Assignment 3:  Start by identifying each section in the code below.  Then annotate as many lines as you can using the 
##comment symbol (#).  If you need help, try googling your question.  Also use the following references:  

###https://www.statmethods.net/r-tutorial/index.html
###https://rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf
###https://rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf
##Don't worry about sections 5 and 6 for now.
##Also - to make the P value comparison you will need to install the package ggpubr.  Let me know if you have issues.


##loading programs####
library(cowplot)
library(ggplot2)
library(tidyverse)
install.packages("ggpubr")
library(ggpubr)

##defining functions####
se <- function(x) sqrt(var(x, na.rm=TRUE)/length(x))


data_summary<-function(x){
  m<-median(x)
  ymin<-m-se(x)
  ymax<-m+se(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

data_summary1<-function(x){
  m<-median(x)
  ymin<-m-(IQR(x)/2)
  ymax<-m+(IQR(x)/2)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

data_summary2<-function(x){
  m<-mean(x)
  ymin<-m-se(x)
  ymax<-m+se(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}


##data summarydataframdefintions####
analyzer<-function(data,gene,control,date,plot_arg,outdir){
  if (dir.exists(outdir)==F) {dir.create(outdir)} 
  df<-read.csv(data, header = TRUE)
  df_1<-df[which(df$censor=="n"),]
  df_2<-df_1[which(df_1$gene_tested==gene),]
  df_3<-df_2[which(df_2$experiment_date==date),]
  df_4c<-df_3[which(df_3$condition==control),]
  df_4e<-df_3[which(df_3$condition==gene),]
  sink(paste0(outdir,"/stat_report.txt"))
  cat(paste0("Stat Report for ",gene," Vs ",control," on ",date,"\n\n"))
  print(df_3 %>% group_by(condition) %>% summarise(n=n(), mean = mean(blinded_score), sd = sd(blinded_score), se = se(blinded_score), median = median(blinded_score), IQR = IQR(blinded_score)))
  con_bs<-df_4c$blinded_score
  exp_bs<-df_4e$blinded_score
  if (length(con_bs)>1 && length(exp_bs>1)) {
    print(ks.test(con_bs,exp_bs))
    print(t_equal<-t.test(con_bs,exp_bs,alternative="two.sided",paired=FALSE,var.equal=TRUE))
    print(t_unequal<-t.test(con_bs,exp_bs,alternative="two.sided",paired=FALSE,var.equal=FALSE))
    print(wilcox<-wilcox.test(con_bs,exp_bs,alternative="two.sided"))
  }
  sink()
  df_3$condition<-relevel(df_3$condition,control)
  plot<-ggplot(data=df_3, aes(x=condition, y=blinded_score, fill = condition))
  plot<-plot+
    geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03)) + 
    scale_fill_manual(values=c("#3C5488","#DC0000")) + 
    scale_color_manual(values=c("black","black"))+
    stat_summary(fun.data=data_summary2, color="black", size=0.5, width=0.3, alpha=0.3,geom="crossbar")+
    stat_compare_means(label = "p.signif", comparisons = list(c(control,gene)))+
    theme_cowplot(font_size = 11)+
    labs(y="Blinded Score", x=NULL)+
    ylim(c(0,1.1*max(df_3$blinded_score)))+
    theme(legend.position = "none")
  if (plot_arg==TRUE){
    save_plot(plot = plot, filename=paste0(outdir,"/",gene,".pdf"),base_width = 3, base_height = 2.7)
  }
  #return(df_3)
}

##plottingdatafordlx4avscontrolon1day####
analyzer(data = "data2.csv", gene = "dlx4a", control = "mcs",date = 20181102, plot_arg = TRUE,outdir = "dlx4a")

##delayre-encodingofstrings###
toc<-read.csv("data2.csv", header = TRUE, stringsAsFactors = F)
toc<-toc %>% select(experiment_date,gene_tested) %>% distinct()

##plottingdataforeachgeneoneachday###
for (i in 1:nrow(toc)) {
  analyzer(data = "data2.csv", gene = toc[i,2], control = "mcs", date = toc[i,1], plot_arg = T, outdir = paste0(toc[i,2],"_",toc[i,1]))
}






