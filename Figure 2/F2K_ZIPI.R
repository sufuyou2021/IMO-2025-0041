setwd("./")
################################################################################
# Load Packages
suppressMessages(library("foreach"))
suppressMessages(library("doMC"))
suppressMessages(library("iterators"))
suppressMessages(library("doParallel"))
suppressMessages(library("parallel"))
suppressMessages(library("Matrix"))
suppressMessages(library("bigmemory"))
suppressMessages(library("biganalytics"))
suppressMessages(library("gRbase"))
suppressMessages(library("gplots"))
suppressMessages(library("ggplot2"))
suppressMessages(library("grid"))
suppressMessages(library("gridExtra"))
suppressMessages(library("data.table"))
suppressMessages(library("plyr"))
suppressMessages(library("ape"))
suppressMessages(library("phyloseq"))
suppressMessages(library("vegan"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("igraph"))
suppressMessages(library("qgraph"))
################################################################################
################################################################################

argv <- commandArgs(TRUE)

  if (length(argv) !=5 ){
    stop('The script needs  files:otu_table; sample_metadata; result_outputfile; data_inputfile ;core_number, please put them in a dir  ')
  }
 
#argv=c("data/16s.uparse_no_pynast_failures_even40000.csv", "data/sample_data_Treat1.tsv", "data/taxonomy.txt", "result", "data")

################################################################################
################################################################################

#Preallocate global data structures 
PARAM <- list();
PARAM$folder.input <- getwd();
PARAM$folder.data <- argv[5]
PARAM$folder.output <- argv[4]              #PUT RESULTS/OUTPUT FOLDER HERE
PARAM$file.otu_data <- argv[3]
PARAM$file.functions <- "/PERSONALBIO/work/microbio/m13/bin/script/functions.community_similarity.R"
PARAM$file.sample_metadata <- argv[2]
PARAM$file.otu_table <- argv[1]
PARAM$cor.use <- "na.or.complete";
PARAM$p.adjust.method <- "hochberg";
PARAM$sample.steps <- c(0.01, 0.02, 0.05, 0.1, 0.15, seq(0.2, 0.9, by=0.1));
PARAM$use.cores <- 1

###################################
commond=paste("mkdir -p ",PARAM$folder.output ,sep="") 
system(commond)
###################################
#Import functions
source(PARAM$file.functions)
###################################
#Minimum OTU size required across all samples
PARAM$thresh.otu_size <- 10;
#Minimum sample size required
PARAM$thresh.sample_size <- 1000;
#Minimum relative number of sequences retained
PARAM$thresh.retained_seqs <- 0.5;
####################################


###################################
#Load data
load(file=paste(PARAM$folder.data, "/All_samples.otu_table.RData", sep=""));
load(file=paste(PARAM$folder.output, "/All_samples.SparCC_global.RData", sep=""));
#Read sample raw metadata
sample.data.raw <- read.table(PARAM$file.sample_metadata, header=T, comment.char = "@", sep="\t", colClasses = "character");
rownames(sample.data.raw) <-as.character(sample.data.raw[,1]);
idx=rownames(sample.data.raw) %in% colnames(ot.2)
sample.data=sample.data.raw[idx,]
samples=rownames(sample.data)
otu.data=read.table(PARAM$file.otu_data, header=T, sep="\t");
rownames(otu.data) <-as.character(otu.data$OTU);
sample.data$Group=factor(as.character(sample.data$Group),levels=unique(sample.data$Group))


############################
############################
#Plot larger network with context
############################
network.samples=samples
#Get current pruned OTU table
curr.ot <- ot.2[which(rowSums(ot.2[,network.samples]) > 0 ), network.samples];
curr.otus <- rownames(curr.ot);

############################
#Generate network plots
############################
#Get SparCC cooc between current OTUs
curr.sparcc <- RMT.sparcc[curr.otus, curr.otus];
curr.sparcc2 <- curr.sparcc
curr.sparcc2[curr.sparcc2 < 0] <- 0;
#Coerce into igraph object
curr.graph <- graph.adjacency(curr.sparcc, mode="upper", weighted=T, diag=F);
curr.graph2 <- graph.adjacency(curr.sparcc2, mode="upper", weighted=T, diag=F);
remove.vertices <- which(degree(curr.graph) == 0);
keep.vertices <- which(degree(curr.graph) >= 1);
curr.graph <- delete.vertices(curr.graph, remove.vertices);
curr.graph2 <- delete.vertices(curr.graph2, remove.vertices);
E(curr.graph)$sim=E(curr.graph)$weight
E(curr.graph)$weight=abs(E(curr.graph)$weight)
curr.ot.rel.raw=data.frame(dgCMatrix2matrix(ot.rel.raw[curr.otus[keep.vertices],samples]))
curr.ot=data.frame(dgCMatrix2matrix(ot[curr.otus[keep.vertices],samples]))

####################################
####################################
#for ggraph 
####################################
####################################
curr.graph.dat <- list();
curr.graph.dat$size <- rowSums(ot[curr.otus[keep.vertices],samples]);

curr.graph.dat$node <- rownames(ot[curr.otus[keep.vertices],samples]);

#get first 9 phylum
temp=aggregate(curr.ot.rel.raw, by=list(otu.data[rownames(curr.ot.rel.raw),"phylum"]), sum)
rownames(temp)=temp$Group.1
temp$Group.1=NULL
taxa10=names(sort(rowSums(temp),decreasing=T))[1:10]

#Color by taxonomy
curr.graph.dat$color.tax <- rep.int("#bdbdbd", length(curr.otus[keep.vertices]));
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[1]] <- "#8DD3C7";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[2]] <- "#FFFFB3";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[3]] <- "#BEBADA";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[4]] <- "#FB8072";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[5]] <- "#80B1D3"
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[6]] <- "#FDB462";
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[7]] <-"#B3DE69"
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[8]] <-"#BC80BD"
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[9]] <-"#CCEBC5"
curr.graph.dat$color.tax[otu.data[curr.otus[keep.vertices], "phylum"] == taxa10[10]] <-"#FFED6F"


col_p=colorRampPalette(brewer.pal(9, "Reds"),bias=0.3,alpha=0.3)(100)
col_n=rev(colorRampPalette(brewer.pal(9, "Greens"),bias=0.3,alpha=0.3)(100))
col_edge=c(col_n,col_p)


############################
############################
#Compute a layout

E(curr.graph)$curved <- 0.25
E(curr.graph)$color <- 0
E(curr.graph)$color=col_edge[100+signif(E(curr.graph)$sim*100,1)]

#community_split using co-occurrence network
com = multilevel.community(curr.graph2, weights = E(curr.graph2)$sim)
subgroup = split(com$names, com$membership)
## subgroup
V(curr.graph)$sg = paste("module",com$membership,sep="_")
modu=names(sort(summary(as.factor(V(curr.graph)$sg),maxsum=10000),decreasing=T))
V(curr.graph)$color="grey"
for (i in 1:10) {
V(curr.graph)$color[which(V(curr.graph)$sg==modu[i])] = c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[c(3,5:7,9,12:13,15,20,19)][i]
}

####
node.id <- c(1:length(V(curr.graph)$name))
names(node.id) <- V(curr.graph)$name

edge <- get.edgelist(curr.graph)
e <- cbind(node.id[edge[,1]],node.id[edge[,2]])
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(curr.graph),area=10*(vcount(curr.graph)^2),repulse.rad=(vcount(curr.graph)^3.1))

############################
############################
# filled by modules
pdf(paste(PARAM$folder.output, "/network_module.pdf",sep=""),width=10,height=10,bg="white")
plot(curr.graph, layout=l, 
    fade=TRUE,
    bg="white",
    trans=TRUE,
    vertex.label=NA,
    vertex.size=log2(curr.graph.dat$size*1000000/length(samples))/6, 
#    vertex.shape=sp,
    vertex.color=V(curr.graph)$color,
    vertex.label.cex=0.5,
    edge.arrow.size=0, 
    edge.color=E(curr.graph)$color, 
    edge.width=E(curr.graph)$weight*4-min(E(curr.graph)$weight)*4+0.001, 
    edge.alpha=0.6,
#    vertex.label.family=NA, 
#    vertex.label.color="black",
#    vertex.label.degree=pi/2,
#    vertex.label.dist=rep(0,length(V(curr.graph)$name)),
#    vertex.pie=pie.values, 
#    vertex.pie.border=list(rep("white"),length(colorg)),
#    vertex.pie.color=list(colorg), 
    vertex.frame.color=NA)
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("topright",title="Modules", c(as.character(modu[1:min(10,length(modu))])), inset=c(-0.05,0),
    title.adj = 0,
    pt.bg=c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[c(3,5:7,9,12:13,15,20,19)][1:min(10,length(modu))],bty = "n",col ="black", pch = 22, cex=1.0))
with(data.frame(V(curr.graph)$sg,V(curr.graph)$color),legend("bottomright",title="Similarity", c("-1.0","-0.9","-0.8","-0.7","0.7","0.8","0.9","1.0"), 
    title.adj = 0,
    col=col_edge[100+c(-0.99,-0.9,-0.8,-0.7,0.7,0.8,0.9,0.99)*100],bty = "n", pch = "--", cex=1.0,pt.cex=2.5,pt.lwd=2))
dev.off();

############################
# filled by Phylum 

############################
# filled by Groups

################################################################
################################################################

############################################
############################################

####################################################
# Within-module degree and participation coefficient
####################################################

names(com$membership) <- com$names
edge <- as.data.frame(edge)
edge$V1 <- as.character(edge$V1)
edge$V2 <- as.character(edge$V2)

edge$V1.m <- com$membership[edge$V1]
edge$V2.m <- com$membership[edge$V2]

# compute degree: nodes ==> modules
degree.mat <- table(data.frame(c(edge$V1,edge$V2),c(edge$V2.m,edge$V1.m)))

# compute Z, P parameters
Z_P.dat <- matrix(0,nrow=vcount(curr.graph), ncol=2)
rownames(Z_P.dat) <- V(curr.graph)$name
colnames(Z_P.dat) <- c("Z", "P")

for(i in rownames(Z_P.dat)){
  Z_P.dat[i, "Z"] <- (degree.mat[i, as.character(com$membership[i])] -
                      mean(degree.mat[, as.character(com$membership[i])][names(com$membership)[com$membership == com$membership[i]]]))/
                      sd(degree.mat[, as.character(com$membership[i])][names(com$membership)[com$membership == com$membership[i]]])
  Z_P.dat[i, "P"] <- 1 - sum((degree.mat[i,]/DEG[i])^2)
}
Z_P.dat.0 <- as.data.frame(Z_P.dat)
Z_P.dat <- Z_P.dat[!is.na(Z_P.dat[,1]),] # remove nodes in modules which have only 2 members / nodes connect to each other within a module
Z_P.dat <- as.data.frame(Z_P.dat)
Z_P.dat$module <- com$membership[rownames(Z_P.dat)]
Z_P.dat$Phylum <- as.character(otu.data[rownames(Z_P.dat), "phylum"])
Z_P.dat$Phylum[!Z_P.dat$Phylum %in% taxa10] <- "Other"
Z_P.dat$Phylum <- factor(Z_P.dat$Phylum, levels=c(taxa10,"Other"))
Z_P.dat$size <- log2(curr.graph.dat$size[rownames(Z_P.dat)]*1000000/length(samples))


#col_roles <- c("Peripherals"="#ffdfba", "Connectors"="#baffc9", "Module hubs"="#bae1ff", "Network hubs"="#e3c6f0")
region.label <- data.frame(x = c(0.28, 0.84, 0.12, 0.88), 
                           y = c(min(Z_P.dat$Z)*1.05, min(Z_P.dat$Z)*1.05, max(max(Z_P.dat$Z)*1.05, 3*0.99),  max(max(Z_P.dat$Z)*1.05, 3*0.99)),
                           label = c("Peripherals", "Connectors", "Module hubs", "Network hubs")
                          )

pdf(paste(PARAM$folder.output, "/network_ZiPi_plot.pdf",sep=""),width=10,height=7.5,bg="white")
ggplot(Z_P.dat,aes(P,Z, colour = Phylum)) + # P threshold = 0.62, Z threshold = 2.5
  #geom_rect(aes(fill = "Peripherals", xmin = -Inf, xmax = 0.62, ymin = -Inf, ymax = 2.5), colour = NA) + 
  #geom_rect(aes(fill = "Connectors", xmin = 0.62, xmax = Inf, ymin = -Inf, ymax = 2.5), colour = NA) + 
  #geom_rect(aes(fill = "Module hubs", xmin = -Inf, xmax = 0.62, ymin = 2.5, ymax = Inf), colour = NA) + 
  #geom_rect(aes(fill = "Network hubs", xmin = 0.62, xmax = Inf, ymin = 2.5, ymax = Inf), colour = NA) + 
  ylim(min(Z_P.dat$Z)*1.05, max(max(Z_P.dat$Z)*1.05, 3)) +
  xlim(0, 1) +
  geom_hline(yintercept = 2.5, linetype=2) +
  geom_vline(xintercept = 0.62, linetype=2) +
  #geom_point(colour = "black" , aes(size=size)) + geom_point(aes(size = size*0.8)) + 
  geom_text(data = region.label, color = "black", size = 5, aes(x = x, y = y, label = label)) +
  geom_point(aes(size = size*0.8, color = Phylum)) +
  xlab("Among-module connectivity(Pi)") + ylab("Within-module connectivity(Zi)") + 
  #scale_fill_manual(name = "Roles", values=col_roles) + scale_size_continuous(range=c(0.5,4), guide="none") + 
  scale_size_continuous(range=c(0.5,4), guide="none") +
  scale_color_manual(values=c(
    "#8DD3C7",
    "#BF8686",
    "#BEBADA",
    "#FB8072",
    "#80B1D3",
    "#FDB462",
    "#B3DE69",
    "#BC80BD",
    "#CCEBC5",
    "#FFED6F",
    "#bdbdbd")) +
  theme_bw() + 
  theme(axis.text.y = element_text(color = "black", size = 15),
        axis.text.x = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 18, face = "bold"),
        axis.title.y = element_text(color = "black", size = 18, face = "bold"),
        legend.text=element_text(face="plain",size=15),legend.title = element_text(size=15,face="bold"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(face = "bold")) +
  guides(colour = guide_legend(override.aes = list(size = 5)))
dev.off()

Z_P_list <- data.frame(OTU_ID=rownames(Z_P.dat), Zi_score=Z_P.dat$Z, Pi_score=Z_P.dat$P)
Z_P.dat.0$Roles = "Peripherals"
Z_P.dat.0$Roles[Z_P.dat.0$Z > 2.5 & Z_P.dat.0$P < 0.62] = "Module hubs"
Z_P.dat.0$Roles[Z_P.dat.0$Z < 2.5 & Z_P.dat.0$P > 0.62] = "Connectors"
Z_P.dat.0$Roles[Z_P.dat.0$Z > 2.5 & Z_P.dat.0$P > 0.62] = "Network hubs"
write.table(Z_P_list, paste(PARAM$folder.output, "/network_ZiPi_list.xls",sep=""),sep="\t",row.names=F,quote = F)

 


#####################################################################################
## for gml
node_topo_list$Domain=as.character(otu.data[names(SIZ), "domain"])
node_topo_list$Phylum=as.character(otu.data[names(SIZ), "phylum"])
node_topo_list$Class=as.character(otu.data[names(SIZ), "class"])
node_topo_list$Order=as.character(otu.data[names(SIZ), "order"])
node_topo_list$Family=as.character(otu.data[names(SIZ), "family"])
node_topo_list$Genus=as.character(otu.data[names(SIZ), "genus"])
node_topo_list$Species=as.character(otu.data[names(SIZ), "species"])
node_topo_list$Top50="FALSE"
node_topo_list[names(head(sort(rowSums(curr.ot),decreasing = T),50)), "Top50"] <- "TRUE"
node_topo_list[is.na(node_topo_list)] <- 0

curr.graph <- delete_vertex_attr(curr.graph, "color") %>% 
              delete_vertex_attr("sg") %>% 
              delete_edge_attr("color") %>% 
              delete_edge_attr("curved")
 
node_topo_list$OTU_ID <- as.character(node_topo_list$OTU_ID)
node_topo_list$Module_group <- as.character(node_topo_list$Module_group)
node_topo_list$Roles <- as.character(node_topo_list$Roles)

for(i in colnames(node_topo_list)){
    curr.graph <- set_vertex_attr(curr.graph, i, value = node_topo_list[,i])
}

write_graph(curr.graph, paste(PARAM$folder.output, "/network_large.gml",sep=""), format = "gml")

warnings()
