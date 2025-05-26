argv <- commandArgs(TRUE)

if (length(argv) != 1){
  stop('The script needs 1 files: 1) lefse_results2.txt' )
}
#argv=c("lefse_results2.txt", "bmap")
suppressMessages(library(ggplot2))
suppressMessages(library(metagenomeSeq))
suppressMessages(library(ggrepel))
suppressMessages(library(grid))
suppressMessages(library(dplyr))
map=read.table("bmap",row.names = 1)
##
main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=8),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))
############################
lefse.data=read.table(argv[1], header=F, sep="\t",row.names=1,check.names = F)
lefse.data$V3=gsub("FALSE","F", lefse.data$V3, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
lefse.data2=lefse.data[which(lefse.data[,3]!=""),]


df=arrange(data.frame(group=as.character(lefse.data2$V3),value=lefse.data2$V4,name=as.character(lefse.data2$V6)),group,value)
df$name=factor(as.character(df$name),levels=c(as.character(df$name)))
df$group=factor(as.character(df$group),levels=unique(map$V4))

#df1=df;df1$name=paste(df$name,"1",sep="");df1$group=paste(df$group,"1",sep="")
#df2=df;df2$name=paste(df$name,"2",sep="");df2$group=paste(df$group,"2",sep="")
#df3=df;df3$name=paste(df$name,"3",sep="");df3$group=paste(df$group,"3",sep="")
#df4=df;df4$name=paste(df$name,"4",sep="");df4$group=paste(df$group,"4",sep="")
#df=rbind(df,df1,df2,df3,df4)
#df=arrange(df,group,value)

if (length(unique(df$group)) > 20){
col=rep(c("magenta4","tomato2","darkolivegreen1","royalblue1","firebrick1","palegoldenrod","grey30",
       "hotpink", "goldenrod4","forestgreen","thistle1","cyan2","maroon4","aquamarine","yellow",
       "grey60","darkred","lightpink","peru","lawngreen","indianred3","mediumseagreen","mediumorchid3",
       "orangered","chocolate2","slateblue","grey80","deepskyblue2","green","gray10"),t=100)[1:length(unique(df$group))] 
}else{
col=c(c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[-4][-9][-9][-11],"#B87333", "#8B00FF", "#CCB38C", "#FF2400")[1:length(unique(df$group))]    #Insufficient color
}
#df=df[-2,]
maxL=1
for (i in 1:length(df$name)){
if(maxL<length(strsplit(as.character(df$name),"")[[i]])) maxL=length(strsplit(as.character(df$name),"")[[i]])
}
pdf(file="lefse_effect_size_rank.pdf",height=trunc(nrow(df)/8+2),width=3.5+(maxL)/20)
p<-ggplot(df, aes(x=name, y=value,fill=group))+
  geom_bar(stat="identity",colour="black",show.legend=TRUE)+coord_flip()+
  theme_bw()+theme(legend.key = element_rect(linetype='blank',fill = 'white'),panel.grid = element_blank())+
  scale_fill_manual(values=col,name="Group",guide = "legend")+guides(fill=guide_legend(title="Group",ncol=1))+xlab("")+ylab("LDA Score (log 10)")+main_theme
p
dev.off()









