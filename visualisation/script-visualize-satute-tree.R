suppressPackageStartupMessages({
  if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
  if(!require("ggtree")) install.packages("ggtree"); library(ggtree)
  if(!require("ggpubr")) install.packages("ggpubr"); library(ggpubr)
  if(!require("ggvenn")) install.packages("ggvenn"); library(ggvenn)
  if(!require("tidytree")) install.packages("tidytree"); library(tidytree)
  if(!require("treeio")) install.packages("treeio"); library(treeio)
  if(!require("phangorn")) install.packages("phangorn"); library(phangorn)
  if(!require("ape"))install.packages("ape");library(ape)
})


library(castor)
library(reshape2)
#library("cowplot")
suppressWarnings(suppressMessages(
  if(!require("BiocManager", quietly = TRUE))
  {
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
    BiocManager::install(c("ggtree","treeio"))
  }))

suppressWarnings(suppressMessages(library("ggtree")))
suppressWarnings(suppressMessages(library("treeio")))


# Input files -----------------------------------
args <- commandArgs()
tree_file <-args[6]
output_file <-args[7]
dir <- args[8]
alpha <-as.numeric(args[9]) # significance level
#option <-args[10] # option could be branch_status, p-value, delta



# ## example test/visual_example
# dir="/home/elgert/Desktop/Cassius/Satute2/test/visual_example"
# option="branch_status"
# tree_file=sprintf("%s/example_resultsRate4.nwk", dir)
# output_file=sprintf("%s/Plot_tree_decision_%s.pdf", dir, option)
  
# # plot tree function
# plot_tree<-function(tree,option,output, alpha=0.05){
#   tree_plot <- ggtree(tree,layout="slanted", aes(col=option))+geom_tippoint() 
#   #+geom_tiplab(align=TRUE, col="black", hjust = -.5, size =6)
#   print(tree_plot)
# }




# Main -------------------
tree2<-read.tree(tree_file)
tree<-read.beast.newick(tree_file)
seq_labels<-unlist(tree["tip.label"])
ttable<-tree%>% as_tibble
ttable[(which(ttable[,1]==ttable[,2])), 4]="R"


#if (option == "branch_status")
#{
  tree_bstatus<-ggtree(tree,layout="slanted", ladderize=FALSE, aes(col=Branch_status)) + geom_tiplab(align=TRUE, colour="black",size=5)+
    xlim(0, max(get_all_distances_to_root(tree2))*1.2)+ 
    theme(legend.position='bottom', plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values=c("#4c9b82", "#E69F00"))
  
   plot_file=sprintf("%s/Plot_tree_%s_decision_branch_status.pdf", dir, output_file)
   pdf(plot_file)
   print(tree_bstatus)
   dev.off()

#} else if (option == "p_value"){
   max_pvalue=as.numeric(max(ttable$p_value[!is.infinite(ttable$p_value)], na.rm=TRUE))
   ttable$p_value[is.na(ttable$p_value)]=0
   ttable$round_p <-round(as.numeric(ttable$p_value),digits = 3)
   tree3<-as.treedata(ttable)
   tree_pvalue<-ggtree(tree3,layout="slanted", ladderize=FALSE, aes(col=as.numeric(p_value))) + geom_tiplab(align=TRUE, colour="black",size=5)+
    xlim(0, max(get_all_distances_to_root(tree2))*1.2)+ 
    theme(legend.position='bottom', plot.title = element_text(hjust = 0.5))+
    labs(colour="p-value")+ geom_text(aes(x=branch, label=round_p))+ 
      scale_colour_gradient2(low="#4c9b82", mid="white", high="#E69F00", midpoint=alpha, limits=c(0,alpha*2),guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))
  
  
  plot_file=sprintf("%s/Plot_tree_%s_decision_p_value_alpha_%s.pdf", dir, output_file, alpha)
  pdf(plot_file)
  print(tree_pvalue)
  dev.off()
#}

  #} else if (option == "delta"){
  max_delta=as.numeric(max(ttable$delta[!is.infinite(ttable$p_value)], na.rm=TRUE))
  ttable$delta[is.na(ttable$delta)]=0
  tree4<-as.treedata(ttable)
  tree_delta<-ggtree(tree4,layout="slanted", ladderize=FALSE, aes(col=as.numeric(delta))) + geom_tiplab(align=TRUE, colour="black",size=5)+
    xlim(0, max(get_all_distances_to_root(tree2))*1.2)+ 
    theme(legend.position='bottom', plot.title = element_text(hjust = 0.5))+
    labs(colour="delta") + geom_text(aes(x=branch, label=delta))
    #+scale_colour_gradient2(low="#4c9b82", mid="white", high="#E69F00", midpoint=mid, limits=c(0,max_delta),guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))
  
  
  plot_file=sprintf("%s/Plot_tree_%s_decision_delta.pdf", dir, output_file)
  pdf(plot_file)
  print(tree_delta)
  dev.off()






