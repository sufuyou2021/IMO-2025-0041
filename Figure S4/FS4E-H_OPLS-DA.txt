#!/PERSONALBIO/Work/Mb/Mb07/.conda/envs/picrust2/bin/Rscript
scal = TRUE
normalize = TRUE
imputation = FALSE 
file = "otu_table_L6.txt"
metadata = "map.txt"
scaling = "pareto"
imput = "mean"

suppressMessages(library(muma))
argv <- commandArgs(TRUE)

if (length(argv) > 0){
    if (length(argv) ==2){
        file<-argv[1]
        metadata <-argv[2]
    }
    if (length(argv) ==3){
        file<-argv[1]
        metadata <-argv[2]
        scaling<-argv[3]
    }
    if (length(argv) ==4){
        file<-argv[1]
        metadata <-argv[2]
        scaling<-argv[3]
        imput <-argv[4]
    }
    if (length(argv)>4 | length(argv)==1 ){
        stop('The script needs 1) otu_table_L6.txt 2) map.txt 3) scalling[pareto(default)|auto|range|vast|median] 4) input[mean(default)|minimum|half.minimum|zero]')
    }
}
#argv=c("l6.txt","map.txt")
data<-read.table(file, sep = "\t", header = TRUE)
map <-read.table(metadata, sep = "\t", header = TRUE, quote = "", colClasses = "character",comment.char = "@")
rownames(data)<-data[,1]
data<-data[,-1]

as.data.frame(cbind(rownames(t(data)),t(t(map))[,4],t(data))) ->df
colnames(df)[1]<-"Samples"
colnames(df)[2]<-"Class"



#############ProcessedTable
comp = df
comp.x = comp[, 3:ncol(comp)]
comp.x = cbind(comp[, 2], comp[, 1], comp.x)

x <- comp.x
x.x <- x[, 3:ncol(x)]

if (!scal) {
    scaling = ""
}
dirout = paste(getwd(), "/Preprocessing_Data_", scaling, 
    "/", sep = "")
dir.create(dirout)
if (imputation) {
    y = x.x
    r = is.na(y)
    for (k in 1:ncol(r)) {
        vec = matrix(r[, k], ncol = 1)
        who.miss.rows = which(apply(vec, 1, function(i) {
            any(i)
        }))
        if (length(who.miss.rows) > nrow(y) * 0.8) {
            warning(paste("The variable -", colnames(y)[k], 
              "- has a number of missing values > 80%, therefore has been eliminated", 
              sep = " "))
            y = y[, -k]
        }
    }
    r = is.na(y)
    who.miss.columns = c()
    for (i in 1:nrow(y)) {
        for (j in 1:ncol(y)) {
            if (r[i, j] == TRUE) {
              if (imput == "mean") {
                v2 = matrix(r[, j], ncol = 1)
                who.miss.rows = which(apply(v2, 1, function(i) {
                  any(i)
                }))
                y[i, j] = mean(y[-who.miss.rows, j])
                print(paste("Imputing missing value of variable -", 
                  colnames(y)[j], "- for the observation -", 
                  rownames(y)[i], "- with", imput, "value", 
                  sep = " "))
              }
              else if (imput == "minimum") {
                v2 = matrix(r[, j], ncol = 1)
                who.miss.rows = which(apply(v2, 1, function(i) {
                  any(i)
                }))
                y[i, j] = min(y[-who.miss.rows, j])
                print(paste("Imputing missing value of variable -", 
                  colnames(y)[j], "- for the observation -", 
                  rownames(y)[i], "- with", imput, "value", 
                  sep = " "))
              }
              else if (imput == "half.minimum") {
                v2 = matrix(r[, j], ncol = 1)
                who.miss.rows = which(apply(v2, 1, function(i) {
                  any(i)
                }))
                y[i, j] = min(y[-who.miss.rows, j])/2
                print(paste("Imputing missing value of variable -", 
                  colnames(y)[j], "- for the observation -", 
                  rownames(y)[i], "- with", imput, "value", 
                  sep = " "))
              }
              else if (imput == "zero") {
                v2 = matrix(r[, j], ncol = 1)
                who.miss.rows = which(apply(v2, 1, function(i) {
                  any(i)
                }))
                y[i, j] = 0
                print(paste("Imputing missing value of variable -", 
                  colnames(y)[j], "- for the observation -", 
                  rownames(y)[i], "- with", imput, "value", 
                  sep = " "))
              }
            }
        }
    }
    pwdi = paste(getwd(), "/Preprocessing_Data_", scaling, 
        "/ImputedMatrix.xls", sep = "")
    write.table(y, pwdi,sep="\t",quote=FALSE)
    x.x = y
}
for (i in 1:nrow(x.x)) {
    for (j in 1:ncol(x.x)) {
        if (as.numeric(x.x[i, j]) <= 0) {
            x.x[i, j] = runif(1, 0, 1e-10)
        }
    }
}
x.x = cbind(comp[, 2], x.x)
write.table(x.x, paste(dirout, "CorrectedTable.xls", sep = ""),quote=FALSE,sep="\t")
pwd.c = paste(getwd(), "/Preprocessing_Data_", scaling, "/CorrectedTable.xls", 
    sep = "")
x <- read.table(pwd.c, sep = "\t", header = TRUE)
x.x <- x[, 2:ncol(x)]
k = matrix(x[, 1], ncol = 1)
if (normalize) {
    x.t <- t(x.x)
    x.s <- matrix(colSums(x.t), nrow = 1)
    uni = matrix(rep(1, nrow(x.t)), ncol = 1)
    area.uni <- uni %*% x.s
    x.areanorm <- x.t/area.uni
    x.areanorm = t(x.areanorm)
    write.table(x.areanorm, paste(dirout, "/ProcessedTable.xls", 
        sep = ""),sep="\t",quote=FALSE)
}else {
    write.table(x.x, paste(dirout, "/ProcessedTable.xls", sep = ""),sep="\t",quote=FALSE)
}
if (scal) {
    if (scaling == "Pareto" | scaling == "pareto" | scaling == 
        "P" | scaling == "p") {
        pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
            "/ProcessedTable.xls", sep = "")
        x <- read.table(pwd.n, sep = "\t", header = TRUE)
        x.x <- x
        x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
        all.sd <- matrix(apply(x.areanorm.tc, 2, sd), nrow = 1)
        uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), 
            ncol = 1)
        all.sdm = uni.exp.all %*% all.sd
        all.sqsd = sqrt(all.sdm)
        all.pareto <- x.areanorm.tc/all.sqsd
        write.table(all.pareto, paste(dirout, "/ProcessedTable.xls", 
            sep = ""),sep="\t",quote=FALSE)
    }
    else if (scaling == "Auto" | scaling == "auto" | scaling == 
        "A" | scaling == "a") {
        pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
            "/ProcessedTable.xls", sep = "")
        x <- read.table(pwd.n, sep = "\t", header = TRUE)
        x.x <- x
        x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
        all.sd <- matrix(apply(x.areanorm.tc, 2, sd), nrow = 1)
        uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), 
            ncol = 1)
        all.sdm = uni.exp.all %*% all.sd
        all.auto <- x.areanorm.tc/all.sdm
        write.table(all.auto, paste(dirout, "/ProcessedTable.xls", 
            sep = ""),sep="\t",quote=FALSE)
    }
    else if (scaling == "Vast" | scaling == "vast" | scaling == 
        "V" | scaling == "v") {
        pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
            "/ProcessedTable.xls", sep = "")
        x <- read.table(pwd.n, sep = "\t", header = TRUE)
        x.x <- x
        x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
        all.sd <- matrix(apply(x.areanorm.tc, 2, sd), nrow = 1)
        uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), 
            ncol = 1)
        all.sdm = uni.exp.all %*% all.sd
        sdm2 = all.sdm^2
        colm = matrix(colMeans(x.x), nrow = 1)
        colm.m = uni.exp.all %*% colm
        num = x.areanorm.tc * colm.m
        vast = num/sdm2
        write.table(vast, paste(dirout, "/ProcessedTable.xls", 
            sep = ""),sep="\t",quote=FALSE)
    }
    else if (scaling == "Range" | scaling == "range" | scaling == 
        "R" | scaling == "r") {
        pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
            "/ProcessedTable.xls", sep = "")
        x <- read.table(pwd.n, sep = "\t", header = TRUE)
        x.x <- x
        x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
        range = c()
        for (i in 1:ncol(x.x)) {
            den = c()
            den = max(x.x[, i]) - min(x.x[, i])
            range = matrix(c(range, den), nrow = 1)
        }
        uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), 
            ncol = 1)
        range.m = uni.exp.all %*% range
        all.range = x.areanorm.tc/range.m
        write.table(all.range, paste(dirout, "/ProcessedTable.xls", 
            sep = ""),sep="\t",quote=FALSE)
    }
    else if (scaling == "Median" | scaling == "median" | 
        scaling == "M" | scaling == "m") {
        pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
            "/ProcessedTable.xls", sep = "")
        x <- read.table(pwd.n, sep = "\t", header = TRUE)
        x.x <- x
        x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
        all.med <- matrix(apply(x.areanorm.tc, 2, median), 
            nrow = 1)
        uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), 
            ncol = 1)
        all.sdm = uni.exp.all %*% all.med
        all.med <- x.areanorm.tc/all.sdm
        write.table(all.med, paste(dirout, "/ProcessedTable.xls", 
            sep = ""),sep="\t",quote=FALSE)
    }
}else {
    pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, 
        "/ProcessedTable.xls", sep = "")
    x <- read.table(pwd.n, sep = "\t", header = TRUE)
    x.x <- x
    x.c = scale(x.x, scale = FALSE)
    write.table(x.c, paste(dirout, "/ProcessedTable.xls", sep = ""),sep="\t",quote=FALSE)
}
pwd.scal = paste(getwd(), "/Preprocessing_Data_", scaling, 
    "/ProcessedTable.xls", sep = "")
x <- read.table(pwd.scal, sep = "\t", header = TRUE)
x.x <- x
pc.all <- prcomp(x.x, center = FALSE, scale = FALSE)
p.v <- matrix(((pc.all$sdev^2)/(sum(pc.all$sdev^2))), ncol = 1)
p.i <- round(p.v * 100, 1)
p.z <- matrix(1, nrow(p.i), 1)
p.f <- cbind(p.i, p.z)
dirout.pca = paste(getwd(), "/PCA_Data_", scaling, "/", sep = "")
dir.create(dirout.pca)
write.table(p.f, paste(dirout.pca, "PCA_P", sep = ""),quote=FALSE,sep="\t")
write.table(pc.all$x, paste(dirout.pca, "PCA_ScoreMatrix.xls", 
    sep = ""),quote=FALSE,sep="\t")
write.table(pc.all$rotation, paste(dirout.pca, "PCA_LoadingsMatrix.xls", 
    sep = ""),quote=FALSE,sep="\t")
pwd.score = paste(getwd(), "/PCA_Data_", scaling, "/", "PCA_ScoreMatrix.xls", 
    sep = "")
Score <- read.table(pwd.score, sep = "\t", header = TRUE)
Score.x <- Score
pwd.load = paste(getwd(), "/PCA_Data_", scaling, "/", "PCA_LoadingsMatrix.xls", 
    sep = "")
Loading <- read.table(pwd.load, sep = "\t", header = TRUE)
Loading.x <- Loading
pwd.pvar = paste(getwd(), "/PCA_Data_", scaling, "/", "PCA_P", 
    sep = "")
Pvar <- read.table(pwd.pvar, sep = "\t", header = TRUE)
Pvar.x <- Pvar
#barplot(Pvar.x[, 1], xlab = "Principal Components", ylab = "Proportion of Variance explained", 
#    main = "Screeplot", ylim = c(0, 100))
#scree = paste(dirout.pca, "Screeplot", scaling, ".pdf", sep = "")
#dev.copy2pdf(file = scree)
col = c(as.character(k))

pairs = c()
if (ncol(Score.x) >= 10) {
    pairs = c(10)
}else {
    pairs = c(ncol(Score.x))
}
K = paste(getwd(), "/Preprocessing_Data_", scaling, "/class.xls", 
    sep = "")
write.table(k, K,quote=FALSE,sep="\t")
x.nn = cbind(k, pc.all$x)
sorted = x.nn[order(x.nn[, 1]), ]
g = c()
for (i in 1:nrow(sorted)) {
    if (any(g == sorted[i, 1])) {
        g = g
    }
    else {
        g = matrix(c(g, sorted[i, 1]), ncol = 1)
    }
}
dirout.g = paste(getwd(), "/Groups", sep = "")
dir.create(dirout.g)
for (i in 1:nrow(g)) {
    vuota <- c()
    fin = matrix(rep(NA, ncol(sorted)), nrow = 1)
    for (j in 1:nrow(sorted)) {
        if (sorted[j, 1] == as.character(g)[i]) {
            vuota <- matrix(sorted[j, ], nrow = 1)
            rownames(vuota) = rownames(sorted)[j]
            fin = rbind(fin, vuota)
        }
    }
    nam = paste("r", as.character(g)[i], sep = ".")
    n = matrix(fin[-1, ], ncol = ncol(sorted))
    n.x = matrix(n[, -1], ncol = ncol(sorted) - 1)
    name = as.matrix(assign(nam, n.x))
    outputfileg = paste("r.", as.character(g)[i], ".xls", sep = "")
    write.table(name, paste(dirout.g, outputfileg, sep = "/"), 
        row.names = FALSE,quote=FALSE,sep="\t")
}



#########OPLSDA
pwd.x = paste(getwd(), "/Preprocessing_Data_", scaling, "/ProcessedTable.xls",sep = "")
x = read.table(pwd.x, header = TRUE, sep = "\t")
x.x = x
pwdK = paste(getwd(), "/Preprocessing_Data_", scaling, "/class.xls", 
             sep = "")
k = read.table(pwdK, header = TRUE, sep = "\t")
x.n = cbind(k, x.x)
sorted = x.n[order(x.n[, 1]), ]
k = matrix(sorted[, 1], ncol = 1)
g = c()
for (i in 1:nrow(sorted)) {
    if (any(g == sorted[i, 1])) {
        g = g
    }
    else {
        g = matrix(c(g, sorted[i, 1]), ncol = 1)
    }
}
Y = matrix(rep(NA, nrow(sorted)), ncol = 1)
for (i in 1:nrow(sorted)) {
    for (l in 1:2) {
        if (sorted[i, 1] == l) {
            Y[i, ] = 0
        }
        else {
            Y[i, ] = 1
        }
    }
}
X = as.matrix(sorted[, -1], ncol = ncol(sorted) - 1)
nf = 1
T = c()
P = c()
C = c()
W = c()
Tortho = c()
Portho = c()
Wortho = c()
Cortho = c()
for (j in 1:nf) {
    w = (t(X) %*% Y) %*% solve(t(Y) %*% Y)
    w1 = t(w) %*% w
    w2 = abs(sqrt(w1))
    w = w %*% solve(w2)
    t = (X %*% w) %*% solve(t(w) %*% w)
    t1 = t(t) %*% t
    c = t(Y) %*% t %*% solve(t1)
    c1 = t(c) %*% c
    u = Y %*% c %*% solve(c1)
    u1 = t(u) %*% u
    u2 = abs(sqrt(u1))
    p = (t(X) %*% t) %*% solve(t1)
    wortho = p - w
    wortho1 = t(wortho) %*% wortho
    wortho2 = abs(sqrt(abs(wortho1)))
    wortho = wortho %*% solve(wortho2)
    tortho = X %*% wortho %*% solve(t(wortho) %*% wortho)
    tortho1 = t(tortho) %*% tortho
    portho = t(X) %*% tortho %*% solve(tortho1)
    cortho = t(Y) %*% tortho %*% solve(tortho1)
    X = X - tortho %*% t(portho)
    T = matrix(c(T, t))
    P = matrix(c(P, p))
    C = matrix(c(C, c))
    W = matrix(c(W, w))
    Tortho = matrix(c(Tortho, tortho))
    Portho = matrix(c(Portho, portho))
    Wortho = matrix(c(Wortho, wortho))
    Cortho = matrix(c(Cortho, cortho))
}
T = matrix(T, ncol = nf)
T = scale(T, scale = FALSE, center = TRUE)
P = matrix(P, ncol = nf)
C = matrix(C, ncol = nf)
W = matrix(W, ncol = nf)
Tortho = matrix(Tortho, ncol = nf)
Portho = matrix(Portho, ncol = nf)
Cortho = matrix(Cortho, ncol = nf)
Wortho = matrix(Wortho, ncol = nf)
Xortho = Tortho %*% t(Portho)
max.pc1 = 1.3 * (max(abs(T[, nf])))
max.pc2 = 1.3 * (max(abs(Tortho[, nf])))
lim = c()
if (max.pc1 > max.pc2) {lim = c(-max.pc1, max.pc1)
}else{lim = c(-max.pc2, max.pc2)}

dirout = paste(getwd(), "/OPLS-DA", scaling, "/", sep = "")
dir.create(dirout)
pwdxdef = paste(dirout, "X_deflated.xls", sep = "")
write.table(X, pwdxdef,quote=FALSE,sep="\t")
pwdT = paste(dirout, "TScore_Matrix.xls", sep = "")
write.table(T, pwdT,quote=FALSE,sep="\t")
pwdTortho = paste(dirout, "TorthoScore_Matrix.xls", sep = "")
write.table(T, pwdTortho,quote=FALSE,sep="\t")
s = as.matrix(sorted[, -1], ncol = ncol(sorted) - 1)
p1 = c()
for (i in 1:ncol(s)) {
    scov = cov(s[, i], T)
    p1 = matrix(c(p1, scov), ncol = 1)
}
pcorr1 = c()
for (i in 1:nrow(p1)) {
    den = apply(T, 2, sd) * sd(s[, i])
    corr1 = p1[i, ]/den
    pcorr1 = matrix(c(pcorr1, corr1), ncol = 1)
}
pwdp1 = paste(dirout, "p1_Matrix.xls", sep = "")
write.table(p1, pwdp1,quote=FALSE,sep="\t")
pwdpcorr1 = paste(dirout, "pcorr1_Matrix.xls", sep = "")
write.table(pcorr1, pwdpcorr1,quote=FALSE,sep="\t")
pc.all <- prcomp(X, center = FALSE, scale = FALSE)
p.f <- summary(pc.all)$importance[2:3,]
write.table(p.f, "OPLS-DA_Proportion.xls",quote=FALSE,sep="\t")
write.table(pc.all$x, "OPLS-DA_ScoreMatrix.xls",quote=FALSE,sep="\t")
write.table(pc.all$rotation,"OPLS-DA_LoadingsMatrix.xls",quote=FALSE,sep="\t")

###############plot

suppressMessages(library(permute))
suppressMessages(library(vegan))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(RColorBrewer))
suppressMessages(library(plyr))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))

Group=factor(map[,4],levels=unique(map[,4]))
if (length(levels(factor(Group)))> 20){
col=rep(c("magenta4","tomato2","darkolivegreen1","royalblue1","firebrick1","palegoldenrod","grey30",           
       "hotpink", "goldenrod4","forestgreen","thistle1","cyan2","maroon4","aquamarine","yellow",
       "grey60","darkred","lightpink","peru","lawngreen","indianred3","mediumseagreen","mediumorchid3",
       "orangered","chocolate2","slateblue","grey80","deepskyblue2","green","gray10"),t=100)[1:length(levels(factor(Group)))] 
}else{
col=c(c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[-4][-9][-9][-11],"#B87333", "#8B00FF", "#CCB38C", "#FF2400")[1:length(levels(factor(Group)))]    #Insufficient color
}
#col=c(c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"))[-4][-9][-9][-11],"#001290", "#2E293F", "#9F838A", "#9F6271", "#607C75", "#207E86", "#FF7063", "#F37EA4", "#BF7100", "#000000")[1:length(levels(factor(Group)))]
#cols=c("#0868AC","#43A2CA","#7BCCC4","#BAE4BC","#F0F9E8")
cols=rev(c(brewer.pal(5,"YlOrRd")))
#colsp=rep(cols,each=20)
shapes=c(15,16,17,18,25)
if (length(levels(factor(Group)))<6) {
    shape=shapes[1:length(levels(factor(Group)))]
}else{
    shape=rep(shapes[2],t=length(levels(factor(Group))))
}
contribution=sort(rowSums(abs(pc.all$rotation[,1:2]*summary(pc.all)$importance[2,1:2])),decreasing = TRUE)
contributors=names(contribution)[which(names(contribution)!="Others")]

plot.data <- data.frame(Axis.1=pc.all$x[,1], Axis.2=pc.all$x[,2],Group=Group)
plot.var_explained <- round(100*summary(pc.all)$importance["Proportion of Variance",1:2], digits=1);
plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
plot.df <- merge(plot.data, plot.centroids, by="Group");
segment=data.frame("seg1"=pc.all$rotation[contributors,1],"seg2"=pc.all$rotation[contributors,2],"seg3"=rep(0,t=length(pc.all$rotation[contributors,1])),"seg4"=rep(0,t=length(pc.all$rotation[contributors,2])),"label"=(rownames(pc.all$rotation[contributors,])))
segment$seg1=segment$seg1*plot.var_explained[1]
segment$seg2=segment$seg2*plot.var_explained[2]
scale=mean(max(plot.data$Axis.1)/max(segment$seg1),min(plot.data$Axis.1)/min(segment$seg1),max(plot.data$Axis.2)/max(segment$seg2),min(plot.data$Axis.2)/min(segment$seg2))

n=length(which(contribution>max(contribution)*0.2))
left <- segment[1:n,][segment[1:n,]$seg1 <0,]
right <- segment[1:n,][segment[1:n,]$seg1>=0,]

num_taxa=min(100,length(segment[,1]))
colsp=colorRampPalette(cols)(num_taxa)

a=max(segment[1:num_taxa,]$seg1)-min(segment[1:num_taxa,]$seg1)
b=max(segment[1:num_taxa,]$seg2)-min(segment[1:num_taxa,]$seg2)

A=max(plot.data$Axis.1)-min(plot.data$Axis.1)
B=max(plot.data$Axis.2)-min(plot.data$Axis.2)

plot.hulls$Group=factor(as.character(plot.hulls$Group),levels=unique(map[,4]))
plot.df$Group=factor(as.character(plot.df$Group),levels=unique(map[,4]))
plot.data$Group=factor(as.character(plot.data$Group),levels=unique(map[,4]))

#plot

sp <- ggplot(data = segment[1:num_taxa,], aes(x=seg1, y=seg2)) +
    geom_text_repel(data=left, aes(x=seg1, y=seg2, label=rownames(left)),
                    color="black",size=3,segment.size=0.1,xlim=c(-Inf,0))+ 
    geom_text_repel(data=right, aes(x=seg1, y=seg2, label=rownames(right)), 
                    color="black", size=3,segment.size=0.1,xlim=c(0.1, Inf))+ 
    geom_point(shape=16,color=colsp,size=2)+
    theme_bw()+
    xlab(paste("PC1 [", plot.var_explained[1], "%]", sep="")) +
    ylab(paste("PC2 [", plot.var_explained[2], "%]", sep="")) +
    theme(axis.text = element_blank(),
          axis.title.x=element_text(color = "black",size = 12),
          axis.title.y=element_text(color = "black", size = 12),
          #axis.ticks.length = unit(.1, "cm"),
          axis.ticks = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
          strip.text = element_text(face = "bold")
    )+ 
    geom_vline(xintercept = 0,linetype=1)+
    geom_hline(yintercept = 0,linetype=1) +
    scale_x_continuous()+ #expand =c(0.05,0.05))+
    scale_y_continuous() #expand =c(0.05,0.05))
    #geom_segment(data= segment[1:5,], aes(xend=seg1*scale, yend=seg2*scale, x=seg3, y=seg4),color="red", size=0.6,alpha=1,arrow =arrow(length = unit(0.2,"cm"))) +
    #geom_text_repel(data=segment[1:20,],aes(x=seg1,y=seg2,label=label),color="black",size=3)+
    #coord_fixed(plot.var_explained[2]/plot.var_explained[1])

Xmin=min(segment[1:num_taxa,]$seg1)
Xmax=max(segment[1:num_taxa,]$seg1)
Ymin=min(segment[1:num_taxa,]$seg2)
Ymax=max(segment[1:num_taxa,]$seg2)
if (((Xmax-Xmin)/(Ymax-Ymin) < 10 && (Xmax-Xmin)/(Ymax-Ymin) > 1) | ((Xmax-Xmin)/(Ymax-Ymin)<1 && (Xmax-Xmin)/(Ymax-Ymin) > 0.1)) {
if (Xmin < Ymin & Xmax < Ymax ){
 sp =sp+ coord_fixed(1,xlim=c(Xmin-abs(Xmax-Ymax)/2,Xmax+abs(Xmax-Ymax)/2),ylim=c(Ymin-abs(Xmin-Ymin)/2,Ymax+abs(Xmin-Ymin)/2),expand =TRUE)
}
if (Xmin > Ymin & Xmax < Ymax ){
  sp =sp+coord_fixed(1,xlim=c(Xmin-abs(Ymax-Ymin-Xmax+Xmin)/2,Xmax+abs(Ymax-Ymin-Xmax+Xmin)/2),ylim=c(Ymin,Ymax),expand =TRUE)
}
if (Xmin < Ymin & Xmax > Ymax ){
  sp =sp+coord_fixed(1,xlim=c(Xmin,Xmax),ylim=c(Ymin-abs(Ymax-Ymin-Xmax+Xmin)/2,Ymax+abs(Ymax-Ymin-Xmax+Xmin)/2),expand =TRUE)
}
if (Xmin > Ymin & Xmax > Ymax ){
  sp =sp+coord_fixed(1,xlim=c(Xmin-abs(Xmin-Ymin)/2,Xmax+abs(Xmin-Ymin)/2),ylim=c(Ymin-abs(Xmax-Ymax)/2,Ymax+abs(Xmax-Ymax)/2),expand =TRUE)
}
}


  
q1 <- ggplot(data = plot.df, aes_string(x="Axis.1", y="Axis.2")) +
    geom_segment(data=plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), color="grey",linetype="solid", size=0.4,alpha=0.9) +
    #geom_vline(xintercept = 0,linetype=2)+
    #geom_hline(yintercept = 0,linetype=2)+
    geom_point(aes(color = Group,fill = Group,shape=Group),size=3,alpha=0.8)+
    #geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.5) +
    theme_bw()+
    xlab(paste("PC1 [", plot.var_explained[1], "%]", sep="")) +
    ylab(paste("PC2 [", plot.var_explained[2], "%]", sep="")) +
    theme(axis.text.x=element_text(color = "black",size = 10),
          axis.text.y=element_text(color = "black", size = 10),
          axis.title.x=element_text(color = "black",size = 12),
          axis.title.y=element_text(color = "black", size = 12),
          axis.ticks.length = unit(.1, "cm"),
          axis.ticks=element_line(color="black"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          #panel.border = element_blank(),
          legend.justification = c(0.5,0.5),
          legend.position="right",
          legend.text=element_text(face="plain",size=12),
          legend.title = element_text(face="plain",size=12),
          legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),
          plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
          strip.text = element_text(face = "bold")
    )+ 
    #geom_text(aes(x=min(plot.data$Axis.1)*3,0,y=0,label=min(plot.data$Axis.1)*3,0),color="black")+
    scale_x_continuous()+ #expand =c(0.05,0.05))+
    scale_y_continuous()+ #expand =c(0.05,0.05))+
    #scale_y_continuous(limits=c(min(plot.data$Axis.2)*3,max(plot.data$Axis.2)*3))+ 
    scale_color_manual(values=col)+
    scale_fill_manual(values=col)+
    scale_shape_manual(values=shape)+
    stat_ellipse(level=0.95, segments=101, alpha=0.8,aes(color=Group),size=0.6,linetype=2) #+
    #geom_segment(data= segment[1:5,], aes(xend=seg1*scale, yend=seg2*scale, x=seg3, y=seg4),color="red", size=0.6,alpha=1,arrow =arrow(length = unit(0.2,"cm"))) +
    #geom_text(data=segment[1:5,],aes(x=seg1*scale*1.1,y=seg2*scale*1.1,label=label),color="black")+
    #coord_fixed(1)

xmin=min(plot.df$Axis.1)
xmax=max(plot.df$Axis.1)
ymin=min(plot.df$Axis.2)
ymax=max(plot.df$Axis.2)
if (((xmax-xmin)/(ymax-ymin) < 10 && (xmax-xmin)/(ymax-ymin) > 1) | ((xmax-xmin)/(ymax-ymin)<1 && (xmax-xmin)/(ymax-ymin) > 0.1)) {
if (xmin < ymin & xmax < ymax ){
 q1 =q1+ coord_fixed(1,xlim=c(xmin-abs(xmax-ymax)/2,xmax+abs(xmax-ymax)/2),ylim=c(ymin-abs(xmin-ymin)/2,ymax+abs(xmin-ymin)/2),expand =TRUE)
}
if (xmin > ymin & xmax < ymax ){
  q1 =q1+coord_fixed(1,xlim=c(xmin-abs(ymax-ymin-xmax+xmin)/2,xmax+abs(ymax-ymin-xmax+xmin)/2),ylim=c(ymin,ymax),expand =TRUE)
}
if (xmin < ymin & xmax > ymax ){
  q1 =q1+coord_fixed(1,xlim=c(xmin,xmax),ylim=c(ymin-abs(ymax-ymin-xmax+xmin)/2,ymax+abs(ymax-ymin-xmax+xmin)/2),expand =TRUE)
}
if (xmin > ymin & xmax > ymax ){
  q1 =q1+coord_fixed(1,xlim=c(xmin-abs(xmin-ymin)/2,xmax+abs(xmin-ymin)/2),ylim=c(ymin-abs(xmax-ymax)/2,ymax+abs(xmax-ymax)/2),expand =TRUE)
}
}
  
p1 <-ggarrange(sp,q1,ncol=2,nrow=1)



q2 <- ggplot(data = plot.df, aes_string(x="Axis.1", y="Axis.2")) +
    geom_segment(data=plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), color="grey",linetype="solid", size=0.4,alpha=0.9) +
    geom_point(aes(color = Group,fill = Group,shape=Group),size=3,alpha=0.8)+
    geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.5) +
    theme_bw()+
    xlab(paste("PC1 [", plot.var_explained[1], "%]", sep="")) +
    ylab(paste("PC2 [", plot.var_explained[2], "%]", sep="")) +
    theme(axis.text.x=element_text(color = "black",size = 10),
          axis.text.y=element_text(color = "black", size = 10),
          axis.title.x=element_text(color = "black",size = 12),
          axis.title.y=element_text(color = "black", size = 12),
          axis.ticks.length = unit(.1, "cm"),
          axis.ticks=element_line(color="black"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.justification = c(0.5,0.5),
          legend.position="right",
          legend.text=element_text(face="plain",size=12),
          legend.title = element_text(face="plain",size=12),
          legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),
          plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
          strip.text = element_text(face = "bold")
    )+ 
    scale_x_continuous()+ #expand =c(0.05,0.05))+
    scale_y_continuous()+ #expand =c(0.05,0.05))+
    #scale_y_continuous(limits=c(min(plot.data$Axis.2)*2,max(plot.data$Axis.2)*2))+
    scale_color_manual(values=col)+
    scale_fill_manual(values=col)+
    scale_shape_manual(values=shape) #+
    #stat_ellipse(level=0.95, segments=101, alpha=0.8,aes(color=Group),size=0.6,linetype=2) +
    #geom_segment(data= segment[1:5,], aes(xend=seg1*scale, yend=seg2*scale, x=seg3, y=seg4),color="red", size=0.6,alpha=1,arrow =arrow(length = unit(0.2,"cm"))) +
    #geom_text(data=segment[1:5,],aes(x=seg1*scale*1.1,y=seg2*scale*1.1,label=label),color="black")+
    #coord_fixed(1)
    
if (((xmax-xmin)/(ymax-ymin) < 10 && (xmax-xmin)/(ymax-ymin) > 1) | ((xmax-xmin)/(ymax-ymin)<1 && (xmax-xmin)/(ymax-ymin) > 0.1)) {
if (xmin < ymin & xmax < ymax ){
 q2 =q2+ coord_fixed(1,xlim=c(xmin-abs(xmax-ymax)/2,xmax+abs(xmax-ymax)/2),ylim=c(ymin-abs(xmin-ymin)/2,ymax+abs(xmin-ymin)/2),expand =TRUE)
}
if (xmin > ymin & xmax < ymax ){
  q2 =q2+coord_fixed(1,xlim=c(xmin-abs(ymax-ymin-xmax+xmin)/2,xmax+abs(ymax-ymin-xmax+xmin)/2),ylim=c(ymin,ymax),expand =TRUE)
}
if (xmin < ymin & xmax > ymax ){
  q2 =q2+coord_fixed(1,xlim=c(xmin,xmax),ylim=c(ymin-abs(ymax-ymin-xmax+xmin)/2,ymax+abs(ymax-ymin-xmax+xmin)/2),expand =TRUE)
}
if (xmin > ymin & xmax > ymax ){
  q2 =q2+coord_fixed(1,xlim=c(xmin-abs(xmin-ymin)/2,xmax+abs(xmin-ymin)/2),ylim=c(ymin-abs(xmax-ymax)/2,ymax+abs(xmax-ymax)/2),expand =TRUE)
}
}

p2 <-ggarrange(sp,q2,ncol=2,nrow=1)


q3 <- ggplot(data = plot.data, aes_string(x="Axis.1", y="Axis.2")) +
    geom_point(aes(color = Group,fill = Group,shape=Group),size=3,alpha=0.8)+
    #geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.5) +
    theme_bw()+
    xlab(paste("PC1 [", plot.var_explained[1], "%]", sep="")) +
    ylab(paste("PC2 [", plot.var_explained[2], "%]", sep="")) +
    theme(axis.text.x=element_text(color = "black",size = 10),
          axis.text.y=element_text(color = "black", size = 10),
          axis.title.x=element_text(color = "black",size = 12),
          axis.title.y=element_text(color = "black", size = 12),
          axis.ticks.length = unit(.1, "cm"),
          axis.ticks=element_line(color="black"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.justification = c(0.5,0.5),
          legend.position="right",
          legend.text=element_text(face="plain",size=12),
          legend.title = element_text(face="plain",size=12),
          legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),
          plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
          strip.text = element_text(face = "bold")
    )+ 
    scale_x_continuous()+ #expand =c(0.05,0.05))+
    scale_y_continuous()+ #expand =c(0.05,0.05))+
    #scale_y_continuous(limits=c(min(plot.data$Axis.2)*2,max(plot.data$Axis.2)*2))+
    scale_color_manual(values=col)+
    scale_fill_manual(values=col)+
    scale_shape_manual(values=shape)+
    #stat_ellipse(level=0.95, segments=101, alpha=0.8,aes(color=Group),size=0.6,linetype=2) +
    #geom_segment(data= segment[1:5,], aes(xend=seg1*scale, yend=seg2*scale, x=seg3, y=seg4),color="red", size=0.6,alpha=1,arrow =arrow(length = unit(0.2,"cm"))) +
    #geom_text(data=segment[1:5,],aes(x=seg1*scale*1.1,y=seg2*scale*1.1,label=label),color="black")+
    geom_vline(xintercept = 0,linetype=2)+
    geom_hline(yintercept = 0,linetype=2) #+
    #coord_fixed(1)

if (((xmax-xmin)/(ymax-ymin) < 10 && (xmax-xmin)/(ymax-ymin) > 1) | ((xmax-xmin)/(ymax-ymin)<1 && (xmax-xmin)/(ymax-ymin) > 0.1)) {
if (xmin < ymin & xmax < ymax ){
 q3 =q3+ coord_fixed(1,xlim=c(xmin-abs(xmax-ymax)/2,xmax+abs(xmax-ymax)/2),ylim=c(ymin-abs(xmin-ymin)/2,ymax+abs(xmin-ymin)/2),expand =TRUE)
}
if (xmin > ymin & xmax < ymax ){
  q3 =q3+coord_fixed(1,xlim=c(xmin-abs(ymax-ymin-xmax+xmin)/2,xmax+abs(ymax-ymin-xmax+xmin)/2),ylim=c(ymin,ymax),expand =TRUE)
}
if (xmin < ymin & xmax > ymax ){
  q3 =q3+coord_fixed(1,xlim=c(xmin,xmax),ylim=c(ymin-abs(ymax-ymin-xmax+xmin)/2,ymax+abs(ymax-ymin-xmax+xmin)/2),expand =TRUE)
}
if (xmin > ymin & xmax > ymax ){
  q3 =q3+coord_fixed(1,xlim=c(xmin-abs(xmin-ymin)/2,xmax+abs(xmin-ymin)/2),ylim=c(ymin-abs(xmax-ymax)/2,ymax+abs(xmax-ymax)/2),expand =TRUE)
}
}


p3 <-ggarrange(sp,q3,ncol=2,nrow=1)

pdf("OPLS-DA.ellipse.pdf",width =15, height =6,bg="white")
p1
dev.off()
pdf("OPLS-DA.hull.pdf",width =15, height =6,bg="white")
p2
dev.off()
pdf("OPLS-DA.common.pdf",width =15, height =6,bg="white")
p3
dev.off()
