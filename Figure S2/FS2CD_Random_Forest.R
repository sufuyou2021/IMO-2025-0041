#!/opt/Software/R/R-3.1.0/bin/Rscript
suppressMessages(library(circlize))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(RColorBrewer))

argv <- commandArgs(TRUE)

if (length(argv) !=3 ){
  stop('The script needs two files:otu_even_depth.xls,importance.tsv, row_number')
}
##argv=c("otu_even_depth.xls","importance.tsv" ,"50")
##argv=c("l5.txt","importance3.tsv","10")
df <- read.table(argv[1], sep = "\t", header = T,row.names=1,check.names = F,comment.char="@")
df2 <- read.table(argv[2], sep = "\t", header = T,check.names = F)
num=as.numeric(argv[3])
df2[,1]=as.character(df2[,1])
df=df[!grepl("unclassified|uncultured|Unknown|unidentified|uncultivated|Other|Incertae", rownames(df)),]
df2=df2[!grepl("unclassified|uncultured|Unknown|unidentified|uncultivated|Other|Incertae", df2[,1]),]
if(nrow(df2)>num){ df2 <- df2[1:num,]}
df1= df[as.character(df2$feature),]

#df1=df
#df1$taxa=factor(rownames(df),levels=c(as.character(rownames(df))))


#mat <- t(df)
mat=df1[,which(colnames(df1)!="taxonomy")]
mat_Z=t(apply(mat, 1, scale))
tax <- colnames(mat_Z)
colnames(mat_Z)=colnames(mat)


if (length(strsplit(rownames(mat_Z)[1],split="")[[1]])> 15 & num > 40){
    df2$ASV_id_used_only_in_this_analysis=1
    for (o in 1:nrow(df2)){
    df2$ASV_id_used_only_in_this_analysis[o]=paste("ASV",o,sep="_")
}
    
    rownames(df2)=df2$feature
    rownames(mat_Z)=df2[rownames(mat_Z),"ASV_id_used_only_in_this_analysis"]
    rownames(df2)=NULL

}
  df2$taxonomy=df1[df2$feature,"taxonomy"]
 write.table(df2,"importance.txt",sep="\t",row.names=F,quote=F)
 



#colnames(mat) <- rep("     ", ncol(mat))
aa=-max(max(mat_Z),-min(mat_Z))
ht = Heatmap(mat_Z,col=colorRamp2(c(-aa,-aa*3/4,-aa*2/4,-aa*1/4, 0, aa*1/4, aa*2/4, aa*3/4, aa), c(c(brewer.pal(5,"YlOrRd")[c(5,4,3,2)]),c(brewer.pal(5,"YlGnBu")[c(1,2,3,4,5)]))),
             show_column_dend = T,width=1.5, cluster_columns=F,cluster_rows=F,row_names_gp = gpar(fontsize = 6),column_names_gp = gpar(fontsize = 8),
             show_heatmap_legend=T,show_row_dend =F,clustering_method_columns = "average",column_dend_reorder = TRUE,
             row_dend_width = unit(4, "cm"),row_names_side="left",column_names_side = "bottom",name = "Z-score")

ha2 = rowAnnotation(Importance = anno_barplot(df2$importance, bar_width = 0.75, gp = gpar(col = "black", fill = "grey"), 
    border =c(Importance= TRUE), which = "row", axis = TRUE),show_annotation_name = T, width = unit(5, "cm"))
#ht + ha2
pdf("feature_importance.pdf",height=3+5*num/50,width=3+(ncol(mat_Z)/2),bg="white")
ht + ha2

dev.off()

