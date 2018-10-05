library(mgcv)

data.files <- commandArgs(TRUE)

cols <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","darkgrey","#A65628")
names(cols) <- c("EastAsia","Africa","WestEurasia","America","Oceania","SouthAsia","CentralAsiaSiberia")


cat("Loading SGDP data...\n")
archaic <- read.table("Archaic_ancestry.txt", as.is=TRUE, header=TRUE)
info <- read.table("Sample_info.txt", as.is=TRUE, header=TRUE)
pca <- read.table("SGDP_PCs.txt", as.is=TRUE, header=TRUE)

cat("Loading PCA weights...\n")
snpweights <- read.table("SGDP_PC_weights.txt.gz", as.is=TRUE, header=TRUE)

cat("Fitting GAM to archiac ancestry...\n")
pop.pca <- aggregate(pca[,1:(NCOL(pca))], mean, by=list(pop=info[rownames(pca),"POP"]))
colnames(pop.pca) <- c("POP", paste0("PC", 1:(NCOL(pop.pca)-1)))
data <- merge(pop.pca, archaic, by="POP")
gam2.n<-gam(Neanderthal~s(PC1,PC2), data=data)
gam4.n<-gam(Neanderthal~s(PC1,PC2,PC3,PC4), data=data)
gam2.d<-gam(Denisovan~s(PC1,PC2), data=data)
gam4.d<-gam(Denisovan~s(PC1,PC2,PC3,PC4), data=data)


cat("Estimating archaic ancestry...\n\n")

n <- length(data.files)
results <- data.frame(Name=rep("",n), Neanderthal=rep(NA, n), Denisovan=rep(NA, n), PC1=rep(NA, n), PC2=rep(NA, n), stringsAsFactors=FALSE)
                              

## load 23andme data
for(i in 1:n){
    data.file <- data.files[i]
    data <- read.table(data.file)
    data$tag <- paste0(data[,2], "_", data[,3])
    data<-data[!duplicated(data$tag),]
    rownames(data) <- data$tag
    name <- strsplit(data.file, "_")[[1]][2]

    cat(paste0(name, "\n"))
    
    data <- data[snpweights[,1],]
    cat(paste0(sum(is.na(data$tag)), " SNPs not in data\n"))
    data <- data[!is.na(data$tag),]
    colnames(data)<-c("ID", "chr", "pos", "bt", "TAG")
    data<-merge(data, snpweights, by="TAG")
    data$gt <- NA
    data$gt <- ifelse(data$bt==mapply(paste0, data$REF, data$REF), 0, data$gt)
    data$gt <- ifelse(data$bt==mapply(paste0, data$REF, data$ALT), 1, data$gt)
    data$gt <- ifelse(data$bt==mapply(paste0, data$ALT, data$REF), 1, data$gt)
    data$gt <- ifelse(data$bt==mapply(paste0, data$ALT, data$ALT), 2, data$gt)
    
    cat(paste0(sum(is.na(data$gt)), " SNPs not typed\n"))
    
    data <- data[!is.na(data$gt),]
    data$gt.rescaled <- (data$gt-data$mu)/sqrt(data$mu * (2-data$mu)/2)
    
    weights <- as.matrix(data[,grepl("^PC", colnames(data))])
    pc<-matrix(data[,"gt.rescaled"],nrow=1) %*% weights / 1e6

    results[i,"Name"] <- name
    results[i,2:5] <- c( round(predict(gam4.n, as.data.frame(pc))*100, 2), round(predict(gam4.d, as.data.frame(pc))*100, 2), pc[1,1], pc[1,2])

    cat(paste0("Neanderthal: ", round(results[i,"Neanderthal"],2), "%\n"))
    cat(paste0("Denisovan: ", round(results[i,"Denisovan"],2), "%\n"))
    
    
}

pdf("Neanderthal_ancestry.pdf")
plot(pca[,1], pca[,2], pch=16, col=cols[info[rownames(pca),"Region"]], xlab=paste0("PC",1), ylab=paste0("PC", 2), main="Neanderthal ancestry")
vis.gam(gam2.n, plot.type="contour",view=c("PC1", "PC2"), add=TRUE, color="bw", nlevels=20)
for(i in 1:n){    
    txt <- paste0( results[i,"Name"], " (", round(results[i,"Neanderthal"],1), "%)")
    points(results[i,"PC1"], results[i,"PC2"], pch=17, cex=1.5, col="black")
    text(results[i,"PC1"], results[i,"PC2"], txt, pos=4)
}
legend("bottomleft", names(cols), col=cols, pch=16, bty="n")
dev.off()

pdf("Denisovan_ancestry.pdf")
plot(pca[,1], pca[,2], pch=16, col=cols[info[rownames(pca),"Region"]], xlab=paste0("PC",1), ylab=paste0("PC", 2), main="Denisovan ancestry")
vis.gam(gam2.d, plot.type="contour",view=c("PC1", "PC2"), add=TRUE, color="bw", nlevels=20)
for(i in 1:n){
    txt <- paste0( results[i,"Name"], " (", round(results[i,"Denisovan"],1), "%)")
    points(results[i,"PC1"], results[i,"PC2"], pch=17, cex=1.5, col="black")
    text(results[i,"PC1"], results[i,"PC2"],txt, pos=4)
}
legend("bottomleft", names(cols), col=cols, pch=16, bty="n")
dev.off()
