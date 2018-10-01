library(UBbipl)
Ocotea.ordered<-data.frame(read.table(file='clipboard', row.names = 1, header=TRUE))
head(Ocotea.ordered)
names(Ocotea.ordered)
attach(Ocotea.ordered)

#1a
VesD<-ordered(VesD)
VesL<-ordered(VesL)
FibL<-ordered(FibL)
RayH<-ordered(RayH)
RayW<-ordered(RayW)
NumVes<-ordered(NumVes)

#1b
datmat <- as.matrix(Ocotea.ordered[,2:7])
Gmat <- indmat(datmat[,1])
for (i in 2:6) Gmat <- cbind(Gmat, indmat(datmat[ ,i]))
Lmat <- diag(apply(Gmat,2,sum))
svd.out <- svd(Gmat %*% solve(sqrt(Lmat))/2)
svd.out$d
round(svd.out$d[-1]^2/sum(svd.out$d[-1]^2),digits=4)

dev.new()
plot(solve(sqrt(Lmat))%*%svd.out$v[,2:3], asp=1,pch=15, col="blue", xlab="Dimension 1", ylab="Dimension 2" )
clp.names <- paste(list("VesD_A","VesD_B","VesD_C","VesD_D","VesD_E","VesL_A","VesL_B","VesL_C","VesL_D","VesL_E","FibL_A","FibL_B","FibL_C","FibL_D","FibL_E","RayH_A","RayH_B","RayH_C","RayH_D","RayH_E",
                        "RayW_A","RayW_B","RayW_C","RayW_D","RayW_E","RayW_F","NumVes_A","NumVes_B","NumVes_C","NumVes_D","NumVes_E","NumVes_F"), sep="")
text(solve(sqrt(Lmat))%*%svd.out$v[,2:3], label=clp.names, pos=3)

#1c
optscores<- solve(sqrt(Lmat))%*%svd.out$v[,2]/sqrt(6)
optscores

#1d
CATPCAbipl(Ocotea.ordered[,2:7], factor.type=c(rep("ord",6)))

Obul <- Ocotea.ordered[,1]=="Obul"
Oken<- Ocotea.ordered[,1]=="Oken"
Opor<- Ocotea.ordered[,1]=="Opor"
col.vec <- rep("green",nrow( Ocotea.ordered))
col.vec[Oken] <- "orange" 
col.vec[Opor]<- "blue"

out <- CATPCAbipl(Ocotea.ordered, plot.samples = 1:nrow(Ocotea.ordered), exp.factor=2, class.vec = Ocotea.ordered[,1],
                  class.cols=c("green", "orange", "blue"), drawbagplots = TRUE, samples.col=col.vec, samples.pch = 8,alpha=95,
                  boxtype = 'o', line.type.bags = rep(2,3), specify.bags = levels(Ocotea.ordered[,1]))

#optimal scores pca
out$zlist

#2a
Ocotea.quantified<-Ocotea.ordered
head(Ocotea.quantified)

Ocotea.quantified$QSpecies[Ocotea.quantified$Species=='Obul']<- 0.09478635
Ocotea.quantified$QSpecies[Ocotea.quantified$Species=='Oken']<- 0.11469250
Ocotea.quantified$QSpecies[Ocotea.quantified$Species=='Opor']<- -0.26985744
Ocotea.quantified$QVesD[Ocotea.quantified$VesD=='A']<--0.28443887
Ocotea.quantified$QVesD[Ocotea.quantified$VesD=='B']<--0.02865714
Ocotea.quantified$QVesD[Ocotea.quantified$VesD=='C']<-0.03400022
Ocotea.quantified$QVesD[Ocotea.quantified$VesD=='D']<-0.16867311
Ocotea.quantified$QVesD[Ocotea.quantified$VesD=='E']<-0.13217469
Ocotea.quantified$QVesL[Ocotea.quantified$VesL=='A']<--0.08669587
Ocotea.quantified$QVesL[Ocotea.quantified$VesL=='B']<--0.28076709
Ocotea.quantified$QVesL[Ocotea.quantified$VesL=='C']<-0.05027416
Ocotea.quantified$QVesL[Ocotea.quantified$VesL=='D']<-0.18821610
Ocotea.quantified$QVesL[Ocotea.quantified$VesL=='E']<-0.12368810 
Ocotea.quantified$QFibL[Ocotea.quantified$FibL=='A']<--0.26993764
Ocotea.quantified$QFibL[Ocotea.quantified$FibL=='B']<--0.05416249
Ocotea.quantified$QFibL[Ocotea.quantified$FibL=='C']<-0.04434368
Ocotea.quantified$QFibL[Ocotea.quantified$FibL=='D']<-0.08574554
Ocotea.quantified$QFibL[Ocotea.quantified$FibL=='E']<-0.20350175
Ocotea.quantified$QRayH[Ocotea.quantified$RayH=='A']<--0.2683576
Ocotea.quantified$QRayH[Ocotea.quantified$RayH=='B']<--0.0762834
Ocotea.quantified$QRayH[Ocotea.quantified$RayH=='C']<-0.0495518
Ocotea.quantified$QRayH[Ocotea.quantified$RayH=='D']<-0.1410263
Ocotea.quantified$QRayH[Ocotea.quantified$RayH=='E']<-0.1683497
Ocotea.quantified$QRayW[Ocotea.quantified$RayW=='A']<--0.25607329
Ocotea.quantified$QRayW[Ocotea.quantified$RayW=='B']<--0.17083084
Ocotea.quantified$QRayW[Ocotea.quantified$RayW=='C']<-0.02679141
Ocotea.quantified$QRayW[Ocotea.quantified$RayW=='D']<-0.15296363
Ocotea.quantified$QRayW[Ocotea.quantified$RayW=='E']<-0.12567075
Ocotea.quantified$QRayW[Ocotea.quantified$RayW=='F']<-0.22240586
Ocotea.quantified$QNumVes[Ocotea.quantified$NumVes=='A']<-0.037181112
Ocotea.quantified$QNumVes[Ocotea.quantified$NumVes=='B']<--0.007348894
Ocotea.quantified$QNumVes[Ocotea.quantified$NumVes=='C']<-0.118516733
Ocotea.quantified$QNumVes[Ocotea.quantified$NumVes=='D']<-0.326090076
Ocotea.quantified$QNumVes[Ocotea.quantified$NumVes=='E']<--0.451713579
Ocotea.quantified$QNumVes[Ocotea.quantified$NumVes=='F']<--0.543686584 

head(Ocotea.quantified)
quant<-Ocotea.quantified[,-(2:7)]
head(quant)
attach(quant)

Obul <- quant[,1]=="Obul"
Oken<- quant[,1]=="Oken"
Opor<- quant[,1]=="Opor"
col.vec <- rep("green",nrow(quant))
col.vec[Oken] <- "orange" 
col.vec[Opor]<- "blue"

pcaout <- PCAbipl(quant[,3:8], G=indmat(quant[,1]),colours=c("green", "orange", "blue"),
                  specify.bags = 1:3, exp.factor = 1.5, marker.size = 0.3, label.size=0.5, pch.samples = rep(20,3))

#classical scaling
AOD.SS(X=quant[,3:8], class.vec=quant[,1], scaled.mat = TRUE,
       X.new.samples=scale( quant[,3:8]), prediction.type="circle", label = FALSE,
       pch.samples.col=c("green","orange","blue"), pch.samples=0:2, pch.means=15:17, dist="Pythagoras", 
       means.plot=TRUE, pch.means.col=c("green","orange","blue"),pch.means.size=1.5,
       expand.markervalsR=rep(1,6),expand.markervalsL=rep(1,6),line.width=1,
       n.int=c(10,30,10,35,20,30),num.points=200,ax=1:6, lwd=2, 
       legend.type=c(TRUE, TRUE, FALSE), legend.x="bottomleft", legend.fractx=0.01, 
       legend.fracty = 0.13)

#biplot with group means

CVAbipl(quant[,3:8], G=indmat(quant[,1]),colours=c("green", "orange", "blue")
        , exp.factor = 1.5, marker.size = 0.3, label.size=0.5, pch.samples = rep(20,3), specify.bags = levels(quant[,1]),
        pch.means = 15:17, colours.means = c("green","orange","blue"), pch.means.size = 1.5)
CVAbipl.pred.regions(quant[,3:8], G=indmat(quant[,1]),colours=c("green", "orange", "blue")
                     , exp.factor = 1.5, marker.size = 0.3, label.size=0.5, pch.samples = rep(20,3), specify.bags = levels(quant[,1]),
                     pch.means = 15:17, colours.means = c("green","orange","blue"), pch.means.size = 1.5,
                     colours.pred.regions = c("lightgreen","yellow","lightblue"), plot.symbol=".", plot.symbol.size = 2.5)

#3a
sites.data<-data.frame(read.table(file='clipboard', header=TRUE))
head(sites.data)
names(sites.data)
attach(sites.data)

#take out missing values
sites.data1<-sites.data[1:11,]
sites.data2<-sites.data[13:37,]
sites.data3<-sites.data[39:129,]
sites.data<-rbind(sites.data1,sites.data2,sites.data3)
sites.data
#metric MDS
sites.euc<-dist(sites.data[,2:12])
cmdscale(sites.euc,k=11,eig=TRUE)
#returns 5 dimensions
max(abs(dist(sites.data[,2:12]) - dist(cmdscale(sites.euc, k = 5))))
#given how large the distance values are, fairly small
cmd<-cmdscale(sites.euc,k=2,eig=TRUE)
plot(cmd$points[,1], cmd$points[,2])

#using optim()
start.config <- cmdscale(sites.euc)
sites.mat <- as.matrix (sites.euc)
y.vec <- as.vector (start.config)

stress.func <- function (y.vec, Delta.mat,f.func=function(x){x}, dim=2)
{ Y <- matrix (y.vec, ncol=dim)
D.mat <- as.matrix (dist(Y))
diff.sq <- (f.func(Delta.mat) - D.mat)^2
diff.sq <- diff.sq[lower.tri(diff.sq)]
sum(diff.sq)/sum(D.mat[lower.tri(D.mat)]^2)  }

optim.out<- optim (y.vec,  stress.func,  Delta.mat=sites.mat)
Y <- matrix(optim.out$par, ncol=2)
optim.out

#plot
plotcol <- sites.data$Site

plotcol <- ifelse(plotcol=="LocA", "black", ifelse(plotcol=="LocB", "blue", ifelse(plotcol=="LocC", "green", ifelse(plotcol=="LocD","yellow","orange"))))

par(pty='s')
plot (Y[,1], Y[,2], asp=1, col=plotcol, pch="+", xlab="Dimension 1", ylab="Dimension 2")

#3b
library(smacof)
S<-smacofSym(sites.euc,type='ordinal')
plot(smacofSym(sites.euc,type='ordinal'),plot.type="confplot", col=plotcol,pch="+", label.conf = list(label="none"))

#3c
#procrustes
fitp<- Procrustes(Y[,1:2], S$conf)
fitp
plot(fitp$Yhat[,1], fitp$Yhat[,2], col=plotcol, xlab="Dimension 1", ylab="Dimension 2")
points(Y[,1],Y[,2], col=plotcol, pch="+")

#4a
climate<-data.frame(read.table(file='clipboard', header=TRUE))
head(climate)
names(climate)
attach(climate)

#missing values
#because not continuous mean didn't make sense
sum(is.na(climate))
climate<-na.omit(climate)

#reverse code
climate$Q8<- 8-climate$Q8
climate$Q9<- 8-climate$Q9
climate$Q10<- 8-climate$Q10
climate$Q11<- 8-climate$Q11
climate$Q12<- 8-climate$Q12
climate$Q13<- 8-climate$Q13
climate$Q14<- 8-climate$Q14
climate$Q15<- 8-climate$Q15
climate$Q16<- 8-climate$Q16
climate$Q17<- 8-climate$Q17
climate$Q18<- 8-climate$Q18
head(climate)

#reliability
library(CTT)
reliability1<-reliability(climate[,2:26])
reliability1
reliability2<-itemAnalysis(climate[,2:26], itemReport=TRUE)
reliability2$itemReport

#sum or mean
climate$score<-rowSums(subset(climate,select=2:26))
climate$score
climate
plot(climate$score)
climate$stscore<-rowSums(subset(climate,select=2:26))/25

ggplot(climate, aes(x=climate$stscore, y=c(0))) + geom_point()

#fa
fa<-factanal(climate[,2:26], factors=1, scores = "Bartlett")
summary(fa)
fa$scores
climate$scores<-fa$scores

library(ggplot2)
ggplot(climate, aes(scores, y=c(0))) + geom_point()
hist(climate$scores)

#5b
#recode data
climate.binary<-climate[,2:26]
attach(climate.binary)
climate.binary[climate.binary==1]<-0
climate.binary[climate.binary==2]<-0
climate.binary[climate.binary==3]<-0
climate.binary[climate.binary==4]<-0
climate.binary[climate.binary==5]<-0
climate.binary[climate.binary==6]<-1
climate.binary[climate.binary==7]<-1
head(climate.binary)
climate.binary<-cbind(climate[,1], climate.binary)

ord<- order(colMeans(climate.binary),decreasing=TRUE)
ord
#best answered items to worst
binary.sorted<-climate.binary[,ord]
binary.sorted
library(ltm)
colMeans(binary.sorted[,2:26])
Tau <- round(-qnorm(colMeans(binary.sorted[,2:26])),2)
Tau
raschfit<-rasch(binary.sorted[,2:26],constraint=cbind(ncol(binary.sorted[,2:26])+1,1))
summary(raschfit)


plot(raschfit, type = "ICC", lwd = 2, legend = TRUE, ncol = 2)
plot(raschfit, type = "IIC", lwd = 2, legend = TRUE, ncol = 2)

climate.binary$raschscores<- (climate.binary$Q7*0.962+climate.binary$Q12*1.002+climate.binary$Q10*1.160+climate.binary$Q18*1.171+climate.binary$Q8*1.181+climate.binary$Q11*1.203+
                                climate.binary$Q4*1.214+ climate.binary$Q24*1.214  +climate.binary$Q1*1.236+climate.binary$Q9*1.236+climate.binary$Q5*1.258+climate.binary$Q15*1.258+
                                climate.binary$Q16*1.258+ climate.binary$Q19*1.258+ climate.binary$Q3*1.269+ climate.binary$Q17*1.281+ climate.binary$Q23*1.281+ climate.binary$Q14*1.292
                              + climate.binary$Q21*1.292+ climate.binary$Q25*1.303+ climate.binary$Q6*1.315+ climate.binary$Q13*1.315+ climate.binary$Q22*1.326+ climate.binary$Q20*1.360
                              + climate.binary$Q2*1.419)
head(climate.binary)

library(ggplot2)
ggplot(climate.binary, aes(raschscores, y=c(0))) + geom_point()
hist(climate.binary$raschscores)

