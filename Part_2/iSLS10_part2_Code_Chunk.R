######################
### Code Chunk 1
######################
getwd() ## verify your current working directory
dir() ## See what files are in the directory
### Confirms that the two main files are in the directory
ldata = read.delim("lipid_data.txt", header=T, as.is=T, check.names=F)
sdata = read.delim("sample_data.txt", header=T, as.is=T, check.names=F)

colnames(ldata) = gsub("\\(-H20\\)", "", colnames(ldata))   
### Remove the neutral loss tag in the lipid names
### For special characters, we need to tell R that they are special by prefixing double backslashes
### Similar commands can be used to manipulate text headers in an automated manner 

ldata[1:5,1:5]
sdata[1:5,]


######################
### Code Chunk 2
######################
### First match the sample IDs between the two data sets
all(ldata$ID %in% sdata$ID)
all(sdata$ID %in% ldata$ID)
all(ldata$ID == sdata$ID) ## All sample IDs are aligned between the two data sets 


######################
### Code Chunk 3
######################
#install.packages("scales")
library(scales)

attach(sdata)
nsamples = nrow(sdata)

hist(BMI, breaks=50)
hist(BMI[DM == 1], breaks=50, add=TRUE, col=2)

plot(SBP ~ Age, pch=19)
plot(SBP ~ Age, cex=.3, pch=19, col=alpha(2, 0.4))  ### red, opaque dots

par(mfrow=c(2,2))
plot(SBP ~ Age, cex=.3, pch=19, col=alpha(2, 0.4))
plot(HbA1c ~ Age, cex=.3, pch=19, col=alpha(2, 0.4))
plot(TG ~ Age, cex=.3, pch=19, col=alpha(2, 0.4))
plot(LDL ~ Age, cex=.3, pch=19, col=alpha(2, 0.4))
dev.off()

## Boxplots
par(mfrow=c(1,2))
boxplot(BMI ~ Gender, boxwex=0.3)
boxplot(HbA1c ~ DM, boxwex=0.3)
dev.off()


######################
### Code Chunk 4
######################
barplot(table(DM), width=0.2, col=c("red","blue"), horiz=TRUE, names.arg = c("Control","Case"))


######################
### Code Chunk 5
######################
## Chi-squared test for categorical data
## and see if there is any correlation / association
tab = table(DM, Gender)
print(tab)
chisq.test(tab)  ### not significant

## Discretizing a continuously scaled variable 
bmi.new = rep(NA, nsamples)
bmi.new[BMI <= 18.5] = 1
bmi.new[BMI > 18.5 & BMI <= 23] = 2
bmi.new[BMI > 23 & BMI <= 27.5] = 3
bmi.new[BMI > 27.5] = 4
tab = table(DM, bmi.new)
print(tab)
chisq.test(tab)  
### significant
detach(sdata)


######################
### Code Chunk 6
######################
rownames(ldata) = ldata$ID ### First column of the data matrix contains identifiers
ldata = ldata[,-1] 
### Since the downstream data analysis will be performed on numerical data only, 
### we put them into ``rownames'' and exclude them from the analysis.

### Check zero or negative concentrations before log transformation
any(ldata <= 0)  
sum(ldata <= 0)  
### Turns out there is a zero or negative value in the data frame somewhere. 

### In case there are zeros, we plug-in a small value for each analyte.
### The downstream analysis should not be sensitive to the choice of this value
### as long as we choose it close to the limite of detection, e.g.
### the smallest positive value across the sample, times 0.9, say. 
### To accomplish this, we run through a "for" loop. 
for(k in 1:ncol(ldata)) {
  ### temporary storage
  xvec = ldata[,k]  
  
  ### take the subset of positive values 
  xpos = xvec[xvec > 0]  
  
  ### get the 90% of the minimum intensity
  ximpute = 0.9 * min(xpos)  
  
  zid = (xvec <= 0)  ### TRUE/FALSE vector in the original vector
  if(any(zid)) {  ### Perform this only if there is a zero
    print(colnames(ldata)[k])  
    xvec[zid] = ximpute  ### Plug in values wherever there is zero
    ldata[,k] = xvec  ### Store back to the original data frame
  }
}

tmp = log2(ldata)


######################
### Code Chunk 7
######################
tmp.pca = prcomp(tmp)
vv = tmp.pca$sdev^2
vv = vv / sum(vv) * 100
vv = round(vv, 2)
print(vv)


######################
### Code Chunk 8
######################
ccc = rep("gray", nsamples)
ccc[sdata$DM == 1] = "red"

Thres = max(abs(tmp.pca$x[,1])) / 2  ### automatically calculate the axis range
XLAB = paste("PC1 (", vv[1], "%)", sep="")  ## PC1
YLAB = paste("PC2 (", vv[2], "%)", sep="")  ## PC2
library(scales)  ### To allow opaqueness in dots

#pdf("PCAplot.pdf", height=5.5, width=5, useDingbats = FALSE)
plot(tmp.pca$x[,1], tmp.pca$x[,2],   ### X=PC1, Y=PC2 coordinate
     col=alpha(ccc,0.2), pch=19,
     xlim=c(-Thres,Thres), ylim=c(-Thres,Thres),   ### Range set equally for both axes
     xlab=XLAB, ylab=YLAB,   ### Automatically assembled above
     main="MEC cohort (N=2,299)", cex=0.8) 
legend("bottomright", c("No event","Incident DM"), pch=19, col=c("gray","red"), cex=0.8)
abline(v=0, lty=2)
abline(h=0, lty=2)
#dev.off()


######################
### Code Chunk 9
######################
############### Draw the entire data (on a relative scale)
tmp.ctr = sweep(tmp, 2, apply(tmp, 2, median))  
### We are normalizing each lipid by its own median value
### so that all lipids are on a comparable scale. 

#install.packages("gplots")
library(gplots)
#pdf("heatmap.pdf", height=20, width=25, useDingbats = FALSE)
heatmap.2(as.matrix(t(tmp.ctr)), trace="n", col=bluered(20), breaks=seq(-2,2,by=0.2), 
          distfun=function(x) as.dist(1-cor(t(x))), 
          hclustfun=function(x) hclust(x, method="average"), 
          ColSideColors = ccc, 
          cexRow=0.2, cexCol=0.2, 
          mar=c(10,10))
#dev.off()


######################
### Code Chunk 10
######################
############### PLS-DA analysis (for fun)
############### Explain why PLS-DA is a "supervised" analysis
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("mixOmics")

library(mixOmics)
sample.class = ifelse(sdata$DM == 1, "Case", "Control")
X = as.matrix(tmp)
Y = sample.class
tmp.out = plsda(X=X, Y=Y, ncomp=2)
#pdf("PLSDA_projection.pdf")
plotIndiv(tmp.out, ind.names = TRUE, ellipse = TRUE, legend = TRUE)
#dev.off()


######################
### Code Chunk 11
######################
head(sdata)


######################
### Code Chunk 12
######################
### We first re-arrange the columns of the sample data matrix
### so that all continuous variables are shifted to the right
sdata = sdata[,c(1:3,5,4,6:14)]  
cp = colnames(sdata)[5:14]  ### Storing variable names

### Now we compute correlations between clinical data and lipids in an one liner
### Recall that ``tmp'' object holds lipidomic data in logarithmic scale (base 2)
cormat = cor(sdata[,5:14], tmp, use="pairwise.complete.obs")

### Plot it in a heatmap
#pdf("cor_heatmap.pdf", height=20, width=6, useDingbats = FALSE)
heatmap.2(as.matrix(t(cormat)), trace="n", main="At baseline",
          col=bluered(20), breaks=seq(-1,1,by=0.1), 
          #distfun=function(x) as.dist(1-cor(t(x))), 
          hclustfun=function(x) hclust(x, method="average"), 
          cexRow=0.2, cexCol=1, 
          mar=c(10,10))
#dev.off()


######################
### Code Chunk 13
######################
### Causal pathway: Lipid --> DM incidence, 
### controlling for HbA1c at baseline, age, gender
### Start the model with no lipid, "base model"
### Then run through a for loop to test contribution of one lipid
sdata$Gender = factor(sdata$Gender, levels=c(1,2))
baseModel = glm(DM ~ Age + Gender + BMI + HbA1c + SBP + HDL + LDL + TG, 
                family=binomial, data=sdata)
summary(baseModel)


######################
### Code Chunk 14
######################
nlipid = ncol(tmp)
nsample = nrow(tmp)
lipid.name = colnames(tmp)
tmp2 = tmp  

for(k in 1:nlipid) {
  mm = mean(tmp[,k], na.rm=TRUE)
  ss = sd(tmp[,k], na.rm=TRUE)
  tmp2[,k] = (tmp[,k] - mm) / ss
}


######################
### Code Chunk 15
######################
### Build model (logistic)
coef.logit = rep(NA, nlipid)
pval.logit = rep(NA, nlipid)

for(k in 1:nlipid) {
  tmpModel = glm(DM ~ Age + Gender + BMI + HbA1c + SBP + HDL + LDL + TG + tmp2[,k], 
                 family=binomial, data=sdata)
  ### Last variable k-th column of tmp2 holds the data for lipid k
  nr = nrow(summary(tmpModel)$coef)
  coef.logit[k] = summary(tmpModel)$coef[nr,1]
  pval.logit[k] = summary(tmpModel)$coef[nr,4]
  if(k %% 50 == 0) print(k)
}
hist(pval.logit, breaks=50)


######################
### Code Chunk 16
######################
# qvalue or Benjamini-Hochberg (BH)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("qvalue")
library(qvalue)
qval.logit = qvalue(pval.logit)$qvalues
plot(pval.logit, qval.logit, cex=.3, pch=19) 
abline(0,1,lty=2) ### properly modeled, it seems

tab.logit = data.frame(lipid=lipid.name, 
                       coefficient=round(coef.logit, 2),
                       pval=pval.logit, qval=qval.logit,
                       stringsAsFactors=FALSE, check.names=FALSE) 
tab.logit = tab.logit[order(tab.logit$pval), ]


######################
### Code Chunk 17
######################
# volcano plot with labels for significant ones
plot(coef.logit, -log10(pval.logit), pch=19, cex=.8, col=alpha("gray", 0.5))
sid = qval.logit <= 0.2
points(coef.logit[sid], -log10(pval.logit[sid]), pch=19, cex=.8, col=alpha("red", 0.5))

tab.logit[tab.logit$qval <= 0.2, ]


######################
### Code Chunk 18
######################
### Step 2: Cox model with time info (time-to-event)
#install.packages("survival")
library(survival)
mid = is.na(sdata$Days)
sdata$Days[mid] = 10000    ### Artificial follow-up date for event-free subjects
baseModelCox = coxph(Surv(Days, DM) ~ Age + Gender + BMI + HbA1c + SBP + HDL + LDL + TG, 
                     data=sdata)
summary(baseModelCox)


######################
### Code Chunk 19
######################
coef.cox = rep(NA, nlipid)
pval.cox = rep(NA, nlipid)

for(k in 1:nlipid) {
  tmpModel = coxph(Surv(Days, DM) ~ Age + Gender + BMI + HbA1c + SBP + HDL + LDL + TG + tmp2[,k], 
                   data=sdata)
  nr = nrow(summary(tmpModel)$coef)
  coef.cox[k] = summary(tmpModel)$coef[nr,1]  ### In Cox, they don't report intercept in summary()
  pval.cox[k] = summary(tmpModel)$coef[nr,5]  ### Also, need to track the columns
  if(k %% 50 == 0) print(k)
}
hist(pval.cox, breaks=100)


######################
### Code Chunk 20
######################
# qvalue or BH
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("qvalue")
library(qvalue)
qval.cox = qvalue(pval.cox)$qvalues
plot(pval.cox, qval.cox, cex=.3, pch=19) 
abline(0,1,lty=2) ### properly modeled, it seems

tab.cox = data.frame(lipid=lipid.name, coefficient=round(coef.cox, 2),
                     pval=pval.cox, qval=qval.cox,
                     stringsAsFactors=FALSE, check.names=FALSE)
tab.cox = tab.cox[order(tab.cox$pval), ]


######################
### Code Chunk 21
######################
# volcano plot with labels for significant ones
plot(coef.cox, -log10(pval.cox), pch=19, cex=.8, col=alpha("gray", 0.5))
sid = (qval.cox <= 0.2)
points(coef.cox[sid], -log10(pval.cox[sid]), pch=19, cex=.8, col=alpha("red", 0.5))

tab.cox[tab.cox$qval <= 0.2, ]


######################
### Code Chunk 22
######################
##### Kaplan-Meier curve
sdata$`SM d16:1/C18:0` = tmp2$`SM d16:1/C18:0`
km1 = survfit(Surv(Days, DM) ~ (`SM d16:1/C18:0` > 0), data=sdata)  ## Fit curve
km1.diff = survdiff(Surv(Days, DM) ~ (`SM d16:1/C18:0` > 0), data=sdata)  ## Test differences

sdata$`Cer d18:0/C18:0` = tmp2$`Cer d18:0/C18:0`
km2 = survfit(Surv(Days, DM) ~ (`Cer d18:0/C18:0` > 0), data=sdata)
km2.diff = survdiff(Surv(Days, DM) ~ (`Cer d18:0/C18:0` > 0), data=sdata)

plot(km1, lty=1, col=c("red","blue"), mark.time=TRUE, 
     xlab="Days", ylab="Survival", main="SM d16:1/C18:0", xlim=c(0,3000))
legend("bottomleft", c("Above zero", "Below zero"), lty=1, col=c("red","blue"))
p.val.1 <- 1 - pchisq(km1.diff$chisq, length(km1.diff$n) - 1)  ## Log-rank test
text(2000, 0.8, paste("p=", round(p.val.1, 4), sep=""))

plot(km2, lty=1, col=c("red","blue"), mark.time=TRUE, 
     xlab="Days", ylab="Survival", main="Cer d18:0/C18:0", xlim=c(0,3000))
legend("bottomleft", c("Above zero", "Below zero"), lty=1, col=c("red","blue"))
p.val.2 <- 1 - pchisq(km2.diff$chisq, length(km2.diff$n) - 1)
text(2000, 0.8, paste("p=", round(p.val.2, 4), sep=""))



######################
### Code Chunk 23
######################
#BiocManager::install("impute")
library(impute)
tmpX = sdata[,5:14]
tmpX.new = impute.knn(as.matrix(tmpX))$data   
### K-nearest neighbor (default K=10)
### Find the most correlated subjects (not variables)
### and estimate the missing data by the average of the 10.

Y = sdata$DM  ### Event 0=noDM / 1=DM
Yt = sdata$Days  ### Follow-up duration or time to event
X = data.frame(tmpX.new, tmp)  ### tmp object holds the lipid data
### In LASSO implementations, each variable is usually standardized within
### But just in case, we'll do that now before running the methods

for(k in 1:ncol(X)) {
  mm = mean(X[,k], na.rm=TRUE)
  ss = sd(X[,k], na.rm=TRUE)
  X[,k] = (X[,k] - mm) / ss
}


######################
### Code Chunk 24
######################
#install.packages("grpreg")
#install.packages("glmnet")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
library(glmnet)
library(grpreg)

### LASSO penalty: impose sparsity constraint on individual lipids
### First run coordinate descent algorithm to get the solution trajectories
### at different penalty levels
logit.lasso = glmnet(x=as.matrix(X), y=Y, family=binomial, alpha=1)
plot(logit.lasso)
### Notice here ``alpha=1'' tells the algorithm to fit LASSO regression
### In LASSO, certain variables are completely suppressed to have zero coefficients
### leading to automatic variable selection
### If alpha is set to 0, the regression is called Ridge regression
### where the coefficients are penalized but not completely obliterated to zero.
### If alpha is set to any value between 0 and 1, it corresponds to Elastic Net (EN) regression. 
### EN is known to have a nice property to deal better with correlated predictors than LASSO.

### Now we must select the optimal penalty level (called lambda) by cross-validation
logit.lasso.cv = cv.glmnet(x=as.matrix(X), y=Y, family=binomial, alpha=1)
plot(logit.lasso.cv)
bestlambda = logit.lasso.cv$lambda.min
bestLASSO = predict(logit.lasso, type="coefficients", s=bestlambda)
bestlambda2 = logit.lasso.cv$lambda.1se  ### 1st dev rule
bestLASSO2 = predict(logit.lasso, type="coefficients", s=bestlambda2)

######################
### Code Chunk 25
######################
### Let's see how E-net is different
logit.enet = glmnet(x=as.matrix(X), y=Y, family=binomial, alpha=0.2)
plot(logit.enet)
logit.enet.cv = cv.glmnet(x=as.matrix(X), y=Y, family=binomial, alpha=0.2)
plot(logit.enet.cv)
bestlambda = logit.enet.cv$lambda.min
bestENET = predict(logit.enet, type="coefficients", s=bestlambda)
bestENET

cbind(bestLASSO, bestENET)


######################
### Code Chunk 26
######################
### Some people may feel uncomfortable with too many lipids of different class of lipids
### being selected. They would like to infuse prior knowledge in the variable selection
var.class = colnames(X)
var.class[grep("Cer.d16.1", colnames(X))] = "Cer_d16_1"
var.class[grep("Cer.d18.0", colnames(X))] = "Cer_d18_0"
var.class[grep("Cer.d18.1", colnames(X))] = "Cer_d18_1"
var.class[grep("Cer.d18.2", colnames(X))] = "Cer_d18_2"
var.class[grep("DHCer.d18.1", colnames(X))] = "DHCer_d18_1"
var.class[grep("GM3.d16.1", colnames(X))] = "GM3_d16_1"
var.class[grep("GM3.d18.1", colnames(X))] = "GM3_d18_1"
var.class[grep("GM3.d18.2", colnames(X))] = "GM3_d18_2"
var.class[grep("MHCer.d16.1", colnames(X))] = "MHCer_d16_1"
var.class[grep("MHCer.d18.0", colnames(X))] = "MHCer_d18_0"
var.class[grep("MHCer.d18.1", colnames(X))] = "MHCer_d18_1"
var.class[grep("MHCer.d18.2", colnames(X))] = "MHCer_d18_2"
var.class[grep("SM.d16.1", colnames(X))] = "SM_d16_1"
var.class[grep("SM.d18.1", colnames(X))] = "SM_d18_1"
var.class[grep("SM.d18.2", colnames(X))] = "SM_d18_2"
var.class[grep("Sphd", colnames(X))] = "Sphd"

levs = unique(var.class)
var.class = factor(var.class, levels=levs)
grpLASSO.fit = grpreg(X=as.matrix(X), y=Y, group=var.class, family="binomial", alpha=1, penalty="grLasso")
plot(grpLASSO.fit)

cv.grpLASSO = cv.grpreg(X=as.matrix(X), y=Y, group=var.class, family="binomial", alpha=1, penalty="grLasso")
plot(cv.grpLASSO)

bestlam.grp = cv.grpLASSO$lambda.min
bestGrpLASSO = coef(grpLASSO.fit, lambda=bestlam.grp)

round( cbind(bestLASSO, bestENET, bestGrpLASSO) , 3)
barplot(bestGrpLASSO[-1], las=2, ylim=c(-0.3,0.3), cex.names=0.3, mar=c(10,10))

######################
### Code Chunk 26
######################

### Random forest
library(randomForest)
set.seed(12345)  ## Fixing the random number generator seed
nvar = ncol(X)
X2 = data.frame(DM=Y, X, stringsAsFactors = FALSE, check.names = FALSE)

X2$DM = factor(X2$DM, levels=c(0,1))  ### Response variable has to be a factor
rf.fit = randomForest(DM ~ ., data=X2, mtry = sqrt(nvar), importance = TRUE, ntree = 1000)
### Obtain the RF classifier
#rf.fit

### See how well it predicts within the training data
yhat.rf = predict(rf.fit, newdata=X2)  ### Making predictions onto the training data itself (no test data available now)
table(yhat.rf, X2$DM)

### Which variables contributed more to the classifier? 
itab = importance(rf.fit)   ### Variable importance measures
ord = order(itab[,3], decreasing=TRUE)
itab = itab[ord, ]
itab[1:20,]  ### Top 20 contributers to the classifier




######################
### Code Chunk 27
######################

### GLASSO for network estimation
#install.packages("huge")
library(huge)
library(gplots)

### Using the package "huge", we will first identify the network of clinical variables and lipids
### To do this, we use a graphical model with sparsity constraint, called Graphical LASSO. 
glasso.out = huge(as.matrix(X), method="glasso")

### There is no equivalent version as "cross-validation" in graphical model estimation. 
### We use a metric called "information criterion". A popular one for GLASSO is RIC
out.select = huge.select(glasso.out, criterion = "ric") ### Rotation information criterion (Lysen 2009)
#plot(out.select)
glasso.fit = out.select$opt.icov  ### taking out the inverse covariance matrix, or precision matrix
rownames(glasso.fit) = colnames(glasso.fit) = colnames(X)

### Export the precision matrix out to a file
glasso.fit = data.frame(Var=colnames(X), glasso.fit, stringsAsFactors=F, check.names=F)
write.table(glasso.fit, "GLASSO.txt", sep="\t", quote=F, row.names=F)

### Visualize the precision matrix
heatmap.2(as.matrix(glasso.fit[,-1]), trace="n", col=bluered(20), breaks=seq(-0.1,0.1,by=0.01))

### Read the precision matrix back in (redundant step)
prec = read.delim("GLASSO.txt", header=T, as.is=T, row.names=1, check.names=F)
### prec refers to precision value

### We now convert precision matrix to partial correlations
### to represent the strength of interaction between pairs of variables
nvar = nrow(prec)  
pcor = prec
### pcor refers to partial correlation

for(i in 1:nvar) {
  for(j in 1:nvar) {
    pcor[i,j] = ifelse(i==j, 1, -1) * prec[i,j] / sqrt(prec[i,i]) / sqrt(prec[j,j])
  }
}

heatmap.2(as.matrix(pcor), trace="n", col=bluered(20), breaks=seq(-0.1,0.1,by=0.01))
### Red represents positive correlation, so now something looks about right. 


###########################################################
#### Create network data: node attributes and edge table
vars = colnames(X)

tvars = c(rep("Clinical",10), rep("Lipids",133))  ### Making class labels for variables
tvars[grep("Cer.", vars)] = "Cer"
tvars[grep("MHCer.", vars)] = "MHCer"
tvars[grep("DHCer.", vars)] = "DHCer"
tvars[grep("GM3.", vars)] = "GM3"
tvars[grep("SM.", vars)] = "SM"
tvars[grep("Sphd", vars)] = "Sphd"

#### Node attribute file
nodes = data.frame(Variable=vars, Type = tvars, stringsAsFactors = FALSE)

#### Now edges and their properties
v1 = v2 = v3 = NULL
for(i in 1:(nvar-1)) {
  for(j in i:nvar) {
    if(pcor[i,j] != 0) {
      v1 = c(v1, vars[i])
      v2 = c(v2, vars[j])
      v3 = c(v3, pcor[i,j])
    }
  }
}
#### Record partial correlations, their absolute values, 
#### and the diretion of correlations (positive / negative)
#### all these will help you visualize the network with more information
edges = data.frame(A=v1, B=v2, pcorr=v3, pcorrAbs=abs(v3), 
                    pcorrSign=ifelse(v3 > 0, "pos", "neg"), 
                    stringsAsFactors = FALSE)
edges = edges[v1 != v2, ]

#### Export the files
write.table(nodes, "Nodes_data.txt", sep="\t", quote=F, row.names=F)
write.table(edges, "Edges_data.txt", sep="\t", quote=F, row.names=F)

#### Off to Cytoscape viz!



