library(ggplot2)
library(phyloseq)
library(dplyr)
library(viridis)
library(cowplot)

library(brms)


PS <- readRDS("tmp/Phyloseq_f.RDS")

## data exploration
PS.ord<-ordinate(PS, "NMDS", "bray")

p1 = plot_ordination(PS, PS.ord, type="taxa", color="Phylum", title="taxa")
print(p1)

head(PS@sam_data)

p2 =plot_ordination(PS, PS.ord, type="samples", color="State")
print(p2)

tax.df <- data.frame(tax_table(PS))
tax.df %>% count(Phylum) -> tax.df


PS.t <- microbiome::transform(PS, "compositional")
ps.df <- psmelt(PS.t)


# ordering the x axis for visualization purposes
PS.t2 <- tax_glom(PS, taxrank="Phylum")

di.m <- vegan::vegdist(PS.t2@otu_table,
                       method="bray")

di.m[is.na(di.m)]<- 0 # defining those as 0 distances
clustering <- hclust(di.m, method="complete")

ps.df$Lab_ID2 <- factor(ps.df$Lab_ID, levels=clustering$labels[clustering$order])
ps.df$Lab_ID <- as.factor(ps.df$Lab_ID)
ps.df$Phylum <- as.factor(ps.df$Phylum)

E1 <- ggplot(tax.df, aes(x=Phylum, y=n, fill=Phylum))+
    geom_col()+
        coord_polar()+
    scale_y_log10("Number of cASVs per phylum")+
    theme_bw()+
    guides(fill="none")

E1



E <- ggplot(ps.df, aes(x=Lab_ID2, y=Abundance, fill=Phylum))+
    geom_bar(stat="identity")+
#    scale_fill_viridis(discrete=T)+
    labs(x="Samples", y="Relative abundance")+
#    facet_grid(.~State, scales = "free_x")+
    theme_classic()+
    guides(fill="none")

E


Fig_abund <- cowplot::plot_grid(E, E1,nrow=1, labels="auto")

ggsave("fig/Exp_composition.pdf", Fig_abund, width=190, height=200, units="mm")

head(PS@sam_data)


############### To do
## HFI is human footprint index
## Infected presence/abscence and load: Infected and Bd.load_co
## deterministic shifts and stochasticity (increase dispersion in diturbed areas and infected animals)
## adults and juveniles
## males and females
## geographical distances, neds to be desintangled with HFI
### predict Bd (classification and regression)

## subset 
PS_w <- subset_samples(PS, State=="Wild")
PS_w

doData <- FALSE

if(doData){
## modeling microbiome composition (b-diversity) by geographical distance, sex, body mass, age,
                                        # infection load and HFI.
metadt <- PS_w@sam_data
key <- sample_names(PS_w)

#1) Jaccard distances
JACM <- as.matrix(phyloseq::distance(PS_w, method="jaccard", type="samples", binary=T))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
# sanity check
all(rownames(JACM)==key)
jac <- c(as.dist(JACM))

#2) Bray distances
BRAM <- as.matrix(phyloseq::distance(PS_w, method="bray", type="samples"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
BRAM <- 1-BRAM
# sanity check
all(rownames(BRAM)==key)
bra <- c(as.dist(BRAM))

## 3) Spatial distance matrix
spa <- c(dist(metadt$Lon, metadt$Lat, method="euclidean"))

## 4) Bodymass distance matrix
bdm <- c(dist(metadt$SMI.rob))

## 5) sex
SexM<-array(0,c(nrow(metadt),nrow(metadt)))
for(i in 1:nrow(metadt)){
    for(j in 1:nrow(metadt)){
        if(metadt$Sex[i]==metadt$Sex[j]){
            SexM[i,j]= "1"
        } else{
            SexM[i,j]= "0"
        }
    }
}
SexM<-c(SexM[lower.tri(SexM)])

## 6) HFI
hfiM <- c(dist(metadt$HFI))

##7) Parasite abundance
infM <- c(dist(metadt$Bd.load_co))

## 8) Age class
#### missing

#Combine these vectors into a data frame
data.dyad<-data.frame(Jac=jac, Bray=bra, spatial=spa, BMI=bdm, Sex=SexM, HFI=hfiM, Infec=infM)

#Now all we need to do is add the identities of both individuals in each dyad as separate columns into the data frame and exclude self-comparisons (as these are not meaningful).
# extracting Individual-combinations present in the matrices
list<-expand.grid(key, key)
# This created individual-to-same-individual pairs as well. Get rid of these:
list<-list[which(list$Var1!=list$Var2),]
# this still has both quantiles in--> add 'unique' key
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(key))
# sanity check that the Individual name combinations are in the same exact order as the lower quantile value vector of the matrices

i=length(key)
JACM[which(rownames(JACM)==list$Var1[i]),which(colnames(JACM)==list$Var2[i])]==jac[i]
# add the names of both individuals participating in each dyad into the data frame

data.dyad$IDA<-list$Var2
data.dyad$IDB<-list$Var1

# Make sure you have got rid of all self comparisons
data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]
#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("spatial","BMI", "HFI", "Infec")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
}

saveRDS(data.dyad, "tmp/data.dyad.RDS")
}

data.dyad <- readRDS("tmp/data.dyad.RDS")

# Now we model Presence/absence
modelJ<-brm(Jac~ 1+ spatial+ BMI + Sex+ HFI+ Infec+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)

saveRDS(modelJ, "tmp/modelJ.rds")

# Now we model Abundances
modelB<-brm(Bray~ 1+ spatial+ BMI + Sex+ HFI+ Infec+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)

saveRDS(modelB, "tmp/modelB.rds")

# ploting, to do

abunM<-summary(modelB)$fixed

res.df<-data.frame(Effect=rownames(abunM),
                   Estimate=abunM$Estimate,
                   lCI=abunM$'l-95% CI',
                   uCI=abunM$'u-95% CI')

res.df

ggplot(res.df, aes(x=Estimate, y=Effect, colour=Effect))+
    geom_point()
    
## random forest regression
library(caret)
library(ranger)
library(pdp)
library(patchwork)
library(microbiome)

PMS.t <- transform(PS_w, "compositional") # transform to make all taxa from 0 to 1
                                        # this is to make it easier to compare among taxa
otu <- PMS.t@otu_table

#PMS.t@tax_table[,6]

colnames(otu) <- paste("ASV", seq(1, length(colnames(otu))), PMS.t@tax_table[,6], sep="_")

df <- data.frame(otu, InfLod=PMS.t@sam_data$Bd.load_co, BMI=PMS.t@sam_data$SMI.rob, Sex=PMS.t@sam_data$Sex, HFI=PMS.t@sam_data$HFI)


set.seed(123)
trainIndex <- createDataPartition(df$InfLod, p=0.8, list=FALSE, times=1)
Inf_df_train <- df[trainIndex,]
Inf_df_test <- df[-trainIndex,]

tgrid <- expand.grid(
    mtry = 1:ncol(otu),
    splitrule = "variance",
    min.node.size = c(5, 10))


set.seed(123)
fitControl <- trainControl(method="repeatedcv", number=10, repeats=10)

#doIgARF <- TRUE

#if (doIgARF){
    
    set.seed(1111)
    rfFit1 <- train(InfLod~., data=Inf_df_train,
                    method="ranger",
                    trControl=fitControl,
                    tuneGrid = tgrid,
                    importance="permutation")

saveRDS(rfFit1, "tmp/rfFit_IgA.rds")

#} else
#    rfFit1 <- readRDS("tmp/rfFit_IgA.rds")

print(rfFit1)
print(rfFit1$finalModel)
test <- predictions <- predict(rfFit1, newdata=IgA <- df <- test)
print(postResample(test <- predictions, IgA <- df <- test$IgA))

pred <- obs <- data.frame(predicted=test <- predictions, observed=IgA <- df <- test$IgA)
cor.test(pred <- obs$predicted, pred <- obs$observed, method="spearman")

corr <-ggplot(pred <- obs, aes(x=predicted, y=observed))+
    geom <- point(size=3, color="black", alpha=0.4)+
                                        #    geom_abline(linetype=5, color="black", size=1, alpha=0.2)+
        stat <- poly <- line(color="black", size=1, alpha=0.2) +
            stat <- poly <- eq() +
                theme <- classic()

imp <- varImp(rfFit1)
imp <- imp$importance
imp$taxa <- rownames(imp)
rownames(imp) <- NULL
imp <- imp[order(-imp$Overall),]
imp20 <- imp[1:20,]
imp20 <- droplevels(imp20)
imp20$taxa <- with(imp20, reorder(taxa, Overall))

IgAt <- data.frame(seq=taxa <- names(PMS)[as.numeric(rownames(imp20))], name=imp20$taxa, importance=imp20$Overall)

library(seqinr)
write.csv2(IgAt, "tmp/IgAtop20.csv")
write.fasta(as.list(IgAt$seq), IgAt$name, "tmp/IgAtop20.fasta")

topImp <-
    ggplot(imp20, aes(y=taxa, x=Overall))+
    geom <- segment( aes(yend=taxa, xend=0)) +
        geom <- point(size=4, color="orange")+
            labs(x="importance", y="")+
            theme <- classic()

Fig4 <- cowplot::plot <- grid(corr, topImp, labels="AUTO", rel <- widths=c(0.6, 1))

top <- features <- imp20$taxa
pd <- plots <- list(NULL)
top <- features <- as.character(top <- features)

for (a in 1:length(top <- features)) {
    pd <- plots[[a]] <-pdp::partial(rfFit1, pred.var=top <- features[a], rug=TRUE)%>%
        autoplot()+
        geom <- hline(yintercept = mean(IgA <- df <- train$IgA), linetype = 2, color = "gray") +
                                        #            scale_y_continuous(limits=c(1.5,2.3)) + # Harmonize the scale of yhat on all plots
            theme(panel.border = element <- rect(colour = "black", fill = NA),
                  panel.background = element <- blank())
    print(paste0("Partial dependence of ", top <- features[a]))
}

fig4 <- 2 <- wrap <- plots(pd <- plots, ncol=4)
fig4 <- cowplot::plot <- grid(Fig4, fig4 <- 2, nrow=2, rel <- heights=c(0.5, 1))library(caret)
library(ranger)
library(pdp)
library(patchwork)

PMS.t <- transform(PMS, "compositional") # transform to make all taxa from 0 to 1
                                        # this is to make it easier to compare among taxa
otu <- PMS.t@otu <- table

colnames(otu) <- paste("ASV", seq(1, length(colnames(otu))), PMS.t@tax <- table[,6], sep="_")
df <- data.frame(otu, IgA=PMS.t@sam <- data$IgA <- inputed)
df$IgA <- log(df$IgA)

set.seed(123)
trainIndex <- createDataPartition(df$IgA, p=0.8, list=FALSE, times=1)
IgA <- df <- train <- df[trainIndex,]
IgA <- df <- test <- df[-trainIndex,]

tgrid <- expand.grid(
    mtry = 1:ncol(otu),
    splitrule = "variance",
    min.node.size = c(5, 10))


set.seed(123)
fitControl <- trainControl(method="repeatedcv", number=10, repeats=10)

doIgARF <- FALSE

if (doIgARF){
    set.seed(1111)
    rfFit1 <- train(IgA~., data=IgA <- df <- train,
                    method="ranger",
                    trControl=fitControl,
                    tuneGrid = tgrid,
                    importance="permutation")
    saveRDS(rfFit1, "tmp/rfFit_IgA.rds")
} else
    rfFit1 <- readRDS("tmp/rfFit_IgA.rds")

print(rfFit1)
print(rfFit1$finalModel)
test <- predictions <- predict(rfFit1, newdata=IgA <- df <- test)
print(postResample(test <- predictions, IgA <- df <- test$IgA))

pred <- obs <- data.frame(predicted=test <- predictions, observed=IgA <- df <- test$IgA)
cor.test(pred <- obs$predicted, pred <- obs$observed, method="spearman")

corr <-ggplot(pred <- obs, aes(x=predicted, y=observed))+
    geom <- point(size=3, color="black", alpha=0.4)+
                                        #    geom_abline(linetype=5, color="black", size=1, alpha=0.2)+
        stat <- poly <- line(color="black", size=1, alpha=0.2) +
            stat <- poly <- eq() +
                theme <- classic()

imp <- varImp(rfFit1)
imp <- imp$importance
imp$taxa <- rownames(imp)
rownames(imp) <- NULL
imp <- imp[order(-imp$Overall),]
imp20 <- imp[1:20,]
imp20 <- droplevels(imp20)
imp20$taxa <- with(imp20, reorder(taxa, Overall))

IgAt <- data.frame(seq=taxa <- names(PMS)[as.numeric(rownames(imp20))], name=imp20$taxa, importance=imp20$Overall)

library(seqinr)
write.csv2(IgAt, "tmp/IgAtop20.csv")
write.fasta(as.list(IgAt$seq), IgAt$name, "tmp/IgAtop20.fasta")

topImp <-
    ggplot(imp20, aes(y=taxa, x=Overall))+
    geom <- segment( aes(yend=taxa, xend=0)) +
        geom <- point(size=4, color="orange")+
            labs(x="importance", y="")+
            theme <- classic()

Fig4 <- cowplot::plot <- grid(corr, topImp, labels="AUTO", rel <- widths=c(0.6, 1))

top <- features <- imp20$taxa
pd <- plots <- list(NULL)
top <- features <- as.character(top <- features)

for (a in 1:length(top <- features)) {
    pd <- plots[[a]] <-pdp::partial(rfFit1, pred.var=top <- features[a], rug=TRUE)%>%
        autoplot()+
        geom <- hline(yintercept = mean(IgA <- df <- train$IgA), linetype = 2, color = "gray") +
                                        #            scale_y_continuous(limits=c(1.5,2.3)) + # Harmonize the scale of yhat on all plots
            theme(panel.border = element <- rect(colour = "black", fill = NA),
                  panel.background = element <- blank())
    print(paste0("Partial dependence of ", top <- features[a]))
}

fig4 <- 2 <- wrap <- plots(pd <- plots, ncol=4)
fig4 <- cowplot::plot <- grid(Fig4, fig4 <- 2, nrow=2, rel <- heights=c(0.5, 1))library(caret)
library(ranger)
library(pdp)
library(patchwork)

PMS.t <- transform(PMS, "compositional") # transform to make all taxa from 0 to 1
                                        # this is to make it easier to compare among taxa
otu <- PMS.t@otu <- table

colnames(otu) <- paste("ASV", seq(1, length(colnames(otu))), PMS.t@tax <- table[,6], sep="_")
df <- data.frame(otu, IgA=PMS.t@sam <- data$IgA <- inputed)
df$IgA <- log(df$IgA)

set.seed(123)
trainIndex <- createDataPartition(df$IgA, p=0.8, list=FALSE, times=1)
IgA <- df <- train <- df[trainIndex,]
IgA <- df <- test <- df[-trainIndex,]

tgrid <- expand.grid(
    mtry = 1:ncol(otu),
    splitrule = "variance",
    min.node.size = c(5, 10))


set.seed(123)
fitControl <- trainControl(method="repeatedcv", number=10, repeats=10)

doIgARF <- FALSE

if (doIgARF){
    set.seed(1111)
    rfFit1 <- train(IgA~., data=IgA <- df <- train,
                    method="ranger",
                    trControl=fitControl,
                    tuneGrid = tgrid,
                    importance="permutation")
    saveRDS(rfFit1, "tmp/rfFit_IgA.rds")
} else
    rfFit1 <- readRDS("tmp/rfFit_IgA.rds")

print(rfFit1)
print(rfFit1$finalModel)
test <- predictions <- predict(rfFit1, newdata=IgA <- df <- test)
print(postResample(test <- predictions, IgA <- df <- test$IgA))

pred <- obs <- data.frame(predicted=test <- predictions, observed=IgA <- df <- test$IgA)
cor.test(pred <- obs$predicted, pred <- obs$observed, method="spearman")

corr <-ggplot(pred <- obs, aes(x=predicted, y=observed))+
    geom <- point(size=3, color="black", alpha=0.4)+
                                        #    geom_abline(linetype=5, color="black", size=1, alpha=0.2)+
        stat <- poly <- line(color="black", size=1, alpha=0.2) +
            stat <- poly <- eq() +
                theme <- classic()

imp <- varImp(rfFit1)
imp <- imp$importance
imp$taxa <- rownames(imp)
rownames(imp) <- NULL
imp <- imp[order(-imp$Overall),]
imp20 <- imp[1:20,]
imp20 <- droplevels(imp20)
imp20$taxa <- with(imp20, reorder(taxa, Overall))

IgAt <- data.frame(seq=taxa <- names(PMS)[as.numeric(rownames(imp20))], name=imp20$taxa, importance=imp20$Overall)

library(seqinr)
write.csv2(IgAt, "tmp/IgAtop20.csv")
write.fasta(as.list(IgAt$seq), IgAt$name, "tmp/IgAtop20.fasta")

topImp <-
    ggplot(imp20, aes(y=taxa, x=Overall))+
    geom <- segment( aes(yend=taxa, xend=0)) +
        geom <- point(size=4, color="orange")+
            labs(x="importance", y="")+
            theme <- classic()

Fig4 <- cowplot::plot <- grid(corr, topImp, labels="AUTO", rel <- widths=c(0.6, 1))

top <- features <- imp20$taxa
pd <- plots <- list(NULL)
top <- features <- as.character(top <- features)

for (a in 1:length(top <- features)) {
    pd <- plots[[a]] <-pdp::partial(rfFit1, pred.var=top <- features[a], rug=TRUE)%>%
        autoplot()+
        geom <- hline(yintercept = mean(IgA <- df <- train$IgA), linetype = 2, color = "gray") +
                                        #            scale_y_continuous(limits=c(1.5,2.3)) + # Harmonize the scale of yhat on all plots
            theme(panel.border = element <- rect(colour = "black", fill = NA),
                  panel.background = element <- blank())
    print(paste0("Partial dependence of ", top <- features[a]))
}

fig4 <- 2 <- wrap <- plots(pd <- plots, ncol=4)
fig4 <- cowplot::plot <- grid(Fig4, fig4 <- 2, nrow=2, rel <- heights=c(0.5, 1))
