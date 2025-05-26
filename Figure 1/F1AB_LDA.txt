suppressMessages(library(phyloseq))
suppressMessages(library(ggtree))
suppressMessages(library(microbiomeViz))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
argv <- commandArgs(TRUE)
if (length(argv) != 4){
      stop('The script needs 3 file and 1 parameter: lefse_table.txt, lefse_tax.txt, lefse_results2.txt and threshold')
}

#argv=c("lefse_table.txt","lefse_tax.txt","lefse_results2.txt","0.0001")   #for test
tax.data=as.matrix(read.table(argv[2], header=T, sep="\t",row.names=1,check.names = F));
otu.data=as.matrix(read.table(argv[1], header=T, sep="\t",row.names=1,comment.char="@",check.names = F));
lefse.data=read.table(argv[3], header=F, sep="\t",row.names=1,check.names = F);
lefse.data$V3=gsub("FALSE","F", lefse.data$V3, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
map=read.table("bmap",row.names = 1)
sample_SUM=colSums(otu.data)
samples=names(sample_SUM)
n=length(sample_SUM)

otu.data=t(t(otu.data)/sample_SUM)
otu.data2=log(otu.data*1000000+1,2)
data=list()
OTU = otu_table(as.matrix(otu.data2), taxa_are_rows = TRUE)
TAX = tax_table(tax.data)
physeq = fix_duplicate_tax(phyloseq(OTU, TAX))
 
tr = parsePhyloseq(physeq, use_abundance = T)
node.data= as.data.frame(ggtree(tr)$data)
rownames(node.data)=node.data$label
taxa=list()
#taxa[[1]]=domain=levels(as.factor(physeq@tax_table[,"domain"]))
taxa[[1]]=phylum=levels(as.factor(physeq@tax_table[,"phylum"]))
taxa[[2]]=class=levels(as.factor(physeq@tax_table[,"class"]))
taxa[[3]]=order=levels(as.factor(physeq@tax_table[,"order"]))
taxa[[4]]=family=levels(as.factor(physeq@tax_table[,"family"]))
taxa[[5]]=genus=levels(as.factor(physeq@tax_table[,"genus"]))
#taxa[[6]]=species=levels(as.factor(physeq@tax_table[,"species"]))
taxa_list=c("phylum","class","order","family","genus")
remain=c();threshold=log(as.numeric(argv[4])*1000000+1,2)
for (j in 5:1){
for ( i in 1:length(taxa[[j]])){
    parent=unique(physeq@tax_table[which(physeq@tax_table[,taxa_list[j]]==taxa[[j]][i]),j])
    member=rownames(physeq@tax_table)[which(physeq@tax_table[,taxa_list[j]]==taxa[[j]][i])]
	if  ( is.na(grep(paste(lefse.data[which(lefse.data[,2]!=""),5],collapse = "|"), physeq@tax_table[member,taxa_list[j]])[1]) && sum(rowSums(physeq@otu_table[member,])/n) < threshold  && (! physeq@tax_table[member,taxa_list[j]] %in% remain ) ) {
	physeq@tax_table[member,taxa_list[j]]=paste("unclassified_",as.character(parent[1]),sep="")
	}else{
	remain=c(remain,parent)
	}
}
}

physeq = fix_duplicate_tax(physeq)
tr = parsePhyloseq(physeq, use_abundance = T)
node.data= as.data.frame(ggtree(tr)$data)
rownames(node.data)=node.data$label
groups=factor(as.character(lefse.data[which(lefse.data[,2]!=""),2]),levels=unique(map$V4)[which( unique(map$V4) %in% as.character(lefse.data[which(lefse.data[,2]!=""),2]))])

if (length(levels(groups)) > 20){
col=rep(c("magenta4","tomato2","darkolivegreen1","royalblue1","firebrick1","palegoldenrod","grey30",
       "hotpink", "goldenrod4","forestgreen","thistle1","cyan2","maroon4","aquamarine","yellow",
       "grey60","darkred","lightpink","peru","lawngreen","indianred3","mediumseagreen","mediumorchid3",
       "orangered","chocolate2","slateblue","grey80","deepskyblue2","green","gray10"),t=100)[1:n] 
}else{
col=c(c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[-4][-9][-9][-11],"#B87333", "#8B00FF", "#CCB38C", "#FF2400")[1:n]    #Insufficient color
}

lefse_lists = data.frame(node=lefse.data[which(lefse.data[,2]!=""),5],
                         color=col[groups],
                         stringsAsFactors = FALSE)
#lefse_lists=lefse_lists[-14,] #bug fix
#lefse_lists=lefse_lists[-15,] #bug fix





p=  ggtree(tr,color="grey50") + theme_tree()+ layout_fan(angle=10)
if (length(lefse.data[which(lefse.data[,2]!=""),5])==0){
pdf(paste("lefse_cladogram.",as.numeric(argv[4]),".pdf",sep=""),width =10, height =10,bg="white")
print(p+geom_point2(aes(size = I(nodeSize)+1), fill = "white", shape = 21,alpha=1))
dev.off()
} else {
short.labs <- c(letters,paste(letters,"1",sep=""),paste(letters,"2",sep=""),paste(letters,"3",sep=""),paste(letters,"4",sep=""),paste(letters,"5",sep=""),paste(letters,"6",sep=""),paste(letters,"7",sep=""),paste(letters,"8",sep=""),paste(letters,"9",sep=""),paste(LETTERS,"1",sep=""),paste(LETTERS,"2",sep=""),paste(LETTERS,"3",sep=""),paste(LETTERS,"4",sep=""),paste(LETTERS,"5",sep=""),paste(LETTERS,"6",sep=""),paste(LETTERS,"7",sep=""),paste(LETTERS,"8",sep=""),paste(LETTERS,"9",sep=""))
get_offset <- function(x) {
    (x * 0.2 + 0.2)^2
}
get_angle <- function(node) {
    data <- p$data
    sp <- tidytree::offspring(data, node)$node
    sp2 <- c(sp, node)
    sp.df <- data[match(sp2, data$node), ]
    mean(range(sp.df$angle))
}
anno.data <- arrange(lefse_lists, node)
hilight.color <- anno.data$color
node_list <- anno.data$node
node_ids <- (p$data %>% filter(label %in% as.character(node_list)) %>% arrange(x))$node
node_list=(p$data %>% filter(label %in% as.character(node_list)) %>% arrange(x))$label
rownames(anno.data)=anno.data$node
anno.data=anno.data[as.character(node_list),]
hilight.color=anno.data[as.character(node_list),"color"]
anno <- rep("white", nrow(p$data))
#######
geom_hilight <- function(node, fill="steelblue", alpha=.5, extend=0, extendto=NULL) {
  
  data = NULL
  stat = "hilight"
  position = "identity"
  show.legend = NA
  na.rm = TRUE
  inherit.aes = FALSE
  check.aes = FALSE
  
  default_aes <- aes_(x=~x, y=~y, node=~node, parent=~parent, branch.length=~branch.length)
  mapping <- default_aes
  
  l <- layer(
    stat=StatHilight,
    data = data,
    mapping = mapping,
    geom = GeomRect,
    position = position,
    show.legend=show.legend,
    inherit.aes = inherit.aes,
    check.aes = check.aes,
    params = list(node=node,
                  fill=fill,
                  alpha=alpha,
                  extend=extend,
                  extendto=extendto,
                  na.rm = na.rm)
  )
  
  return(l)
}
#########
for (i in length(node_ids):1) {
    n <- node_ids[i]
    color <- hilight.color[i]
    anno[n] <- color
    mapping <- p$data %>% filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    offset <- get_offset(nodeClass)
    p <- p + geom_hilight(node = n, fill = color,alpha = 0.4, extendto = 6+0.25*as.numeric(nodeClass))
}

p$layers <- rev(p$layers)
p <- p + geom_point2(aes(size = I(nodeSize)+1), fill = anno, shape = 21,alpha=1)
short.labs.anno <- NULL
anno.depth=7
 for (i in 1:length(node_ids)) {
         n <- node_ids[i]
         mapping <- p$data %>% filter(node == n)
        # nodeClass <- as.numeric(mapping$nodeClass)
         lab <- short.labs[i]
        # short.labs <- short.labs[-1]
         if (is.null(short.labs.anno)) {
             short.labs.anno = data.frame(lab = lab, annot = mapping$label,stringsAsFactors = F)
         }else {
             short.labs.anno = rbind(short.labs.anno, c(lab,mapping$label))
         }
        # offset <- get_offset(nodeClass) - 0.4
        # angle <- get_angle(n) + 90
       # p <- p + geom_cladelabel(node = n, label = lab,  angle = angle, fontsize = 1.5 + sqrt(nodeClass), offset = offset, barsize = 0, hjust = 0.5)
        p$data$label[n]=lab
       
 }
p=p+geom_text2(aes(subset=(node %in% node_ids),label=label,x=x+0.4,angle=ifelse((angle>80 & angle<260) ,angle-190,angle)))

set_hilight_legend <- function(p, color, label, alpha=1) {
	d <- data.frame(color=color, clade=label, x=0, y=1, alpha=alpha)
	p + geom_rect(aes_(fill=~clade, xmin=~x, xmax=~x, ymin=~y, ymax=~y), data=d, inherit.aes=F) +
		guides(fill=guide_legend(title="Taxa",override.aes=list(fill=alpha(d$color, d$alpha)),nrow=25))
}
q = set_hilight_legend(p, as.character(anno.data$color), factor(paste(short.labs.anno$lab,": ",short.labs.anno$annot,sep=""),levels=paste(short.labs.anno$lab,": ",short.labs.anno$annot,sep="")),alpha=1) + theme(legend.position="right",legend.title=element_text(size=18))

pdf(paste("lefse_cladogram.",as.numeric(argv[4]),".pdf",sep=""),width =13+(length(node_ids)%/%25)*3.2, height =10,bg="white", onefile=F)

plot.new()
print(q)
palette(col)
with(data.frame(group=levels(groups),col[1:length(levels(groups))]),legend("topleft",title="Groups",levels(groups), inset=c(0.01,0.1),
    title.adj = 0, pt.bg=c(as.character(col[1:length(levels(groups))])),bty = "n", pch = 22, cex=1.5,col=c(as.character(col[1:length(levels(groups))]))))
dev.off()
}
