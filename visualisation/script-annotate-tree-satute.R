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

library(stringr)
library(ggtree)

# Input files -----------------------------------
args <- commandArgs()
tree_file <-args[6]
res_file <-args[7]
dir <- args[8]
output_file <-args[9]


## example octo-kraken-msa-test
dir="/home/elgert/Desktop/Cassius/Satute2/test/octo-kraken-msa-test/"
tree_file=sprintf("%s/example.phy.treefile", dir)
res_file=sprintf("%s/tree_resultsRate1_table",dir)
output_file=sprintf("%s/outfile",dir)

# Main -------------------

tree<-read.tree(tree_file)
seq_labels<-unlist(tree["tip.label"])

# tree_plot <- ggtree(tree,layout="slanted")+
#   geom_tiplab(align=TRUE, col="black", hjust = -.5, size =6)+geom_tippoint() 
# print(tree_plot)
ttable<-tree%>% as_tibble
ttable[(which(ttable[,1]==ttable[,2])), 4]="R"


res_satute=read.csv(res_file, header=FALSE, sep="\t")
colnames(res_satute)=c("Order", "delta", "c_s", "p_value", "Branch status","T2T status", "Branch")
pairs<-strsplit(str_sub(res_satute$Branch, start=1, end=-2), "*-")
pairs<-lapply(pairs, unlist)
pair_1<-unlist(lapply(pairs,function(x) x[1]))
pair_2<-unlist(lapply(pairs,function(x) x[2]))
pair_1<-gsub(" ", "", pair_1, fixed = TRUE)
pair_2<-gsub(" ", "", pair_2, fixed = TRUE)
res_satute$node_1<-pair_1
res_satute$node_2<-pair_2

info_satute <- res_satute[, c("Branch status","p_value", "node_2")]
colnames(info_satute)<-c("Branch_status","p_value", "label")

table2<-ttable%>% left_join(info_satute, by='label')
tree2<-as.treedata(table2);


#save nexus file 
tree_nexus=sprintf("%s/%s", dir, output_file)
write.beast(tree2, file=tree_nexus)
# tree_nwk=sprintf("%s/%s_treefile_extended.nwk", dir, name)
# write.beast.newick(tree2,tree_nwk)

# plot of tree 
tree_plot <- ggtree(tree2,layout="slanted", aes(col=Branch_status))+geom_tippoint() 
  #+geom_tiplab(align=TRUE, col="black", hjust = -.5, size =6)
print(tree_plot)

plot_file=sprintf("%s/Tree_Rplot.pdf", dir)
pdf(plot_file)
print(tree_plot)
dev.off()



