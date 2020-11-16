### FB 26/03/2019 -- Plots of the variances covariances matrices

rm(list=ls())
graphics.off()

### Loads data

FIM1<-read.csv("../frequentist/FIM1.csv")
FIM2<-read.csv("../frequentist/FIM2.csv")
FIM1<-FIM1[,-1]
FIM2<-FIM2[,-1]

### Construct correlation matrices
library(corrplot)

Sigma = solve(FIM1)
Rho = cov2cor(Sigma) #correlation matrix

## see https://stackoverflow.com/questions/51518618/using-expression-to-label-corrplot for labelling in corrplot
par(mfrow=c(1,1))
rownames(Rho) = c("r",":gamma",":sigma[1]","s","Q",":sigma[2]","C","D",":sigma[3]")
colnames(Rho) = c("r",":gamma",":sigma[1]","s","Q",":sigma[2]","C","D",":sigma[3]")
corrplot(Rho, method="circle")
pdf("InverseFIM1_corrplot.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()

Sigma = solve(FIM2)
Rho = cov2cor(Sigma) #correlation matrix
par(mfrow=c(1,1))
rownames(Rho) = c("r",":gamma",":sigma[1]","s","Q",":sigma[2]","C","D")
colnames(Rho) = c("r",":gamma",":sigma[1]","s","Q",":sigma[2]","C","D")
corrplot(Rho, method="circle")
pdf("InverseFIM2_corrplot.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()

## LC data

FIM1<-read.csv("../frequentist/FIM1_noisyLC.csv")
FIM2<-read.csv("../frequentist/FIM2_noisyLC.csv")
FIM1<-FIM1[,-1]
FIM2<-FIM2[,-1]


Sigma = solve(FIM1)
Rho = cov2cor(Sigma) #correlation matrix
par(mfrow=c(1,1))
rownames(Rho) = c("r",":gamma",":sigma[1]","s","Q",":sigma[2]","C","D",":sigma[3]")
colnames(Rho) = c("r",":gamma",":sigma[1]","s","Q",":sigma[2]","C","D",":sigma[3]")
corrplot(Rho, method="circle")
pdf("InverseFIM1_noisyLC_corrplot.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()

Sigma = solve(FIM2)
Rho = cov2cor(Sigma) #correlation matrix
par(mfrow=c(1,1))
rownames(Rho) = c("r",":gamma",":sigma[1]","s","Q",":sigma[2]","C","D")
colnames(Rho) = c("r",":gamma",":sigma[1]","s","Q",":sigma[2]","C","D")
corrplot(Rho, method="circle")
pdf("InverseFIM2_noisyLC_corrplot.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()


### Reparameterized model


### Loads data

FIM1<-read.csv("../frequentist/FIM1_reparam.csv")
FIM2<-read.csv("../frequentist/FIM2_reparam.csv")
FIM1<-FIM1[,-1]
FIM2<-FIM2[,-1]

### Construct correlation matrices

Sigma = solve(FIM1)
Rho = cov2cor(Sigma) #correlation matrix


par(mfrow=c(1,1))
rownames(Rho) = c("r","K",":sigma[1]","s","Q",":sigma[2]","a","h",":sigma[3]")
colnames(Rho) = c("r","K",":sigma[1]","s","Q",":sigma[2]","a","h",":sigma[3]")
corrplot(Rho, method="circle")
pdf("InverseFIM1_corrplot_reparam.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()

Sigma = solve(FIM2)
Rho = cov2cor(Sigma) #correlation matrix
par(mfrow=c(1,1))
rownames(Rho) = c("r","K",":sigma[1]","s","Q",":sigma[2]","a","h")
colnames(Rho) = c("r","K",":sigma[1]","s","Q",":sigma[2]","a","h")
corrplot(Rho, method="circle")
pdf("InverseFIM2_corrplot_reparam.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()

## LC data

FIM1<-read.csv("../frequentist/FIM1_noisyLC_reparam.csv")
FIM2<-read.csv("../frequentist/FIM2_noisyLC_reparam.csv")
FIM1<-FIM1[,-1]
FIM2<-FIM2[,-1]


Sigma = solve(FIM1)
Rho = cov2cor(Sigma) #correlation matrix
par(mfrow=c(1,1))
rownames(Rho) = c("r","K",":sigma[1]","s","Q",":sigma[2]","a","h",":sigma[3]")
colnames(Rho) = c("r","K",":sigma[1]","s","Q",":sigma[2]","a","h",":sigma[3]")
corrplot(Rho, method="circle")
pdf("InverseFIM1_noisyLC_corrplot_reparam.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()

Sigma = solve(FIM2)
Rho = cov2cor(Sigma) #correlation matrix
par(mfrow=c(1,1))
rownames(Rho) = c("r","K",":sigma[1]","s","Q",":sigma[2]","a","h")
colnames(Rho) = c("r","K",":sigma[1]","s","Q",":sigma[2]","a","h")
corrplot(Rho, method="circle")
pdf("InverseFIM2_noisyLC_corrplot_reparam.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()

