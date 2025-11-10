library(ggplot2)
library(phyloseq)
library(dplyr)
library(viridis)
library(cowplot)
library(caret)
library(ranger)
library(pdp)
library(patchwork)
library(microbiome)
library(brms)
library(ggpmisc)
library(gridExtra)
library(ggmcmc)
library(ggridges)
library(paletteer)
library(viridis)

PS <- readRDS("tmp/Phyloseq_f.RDS")

## description of the environmental variables
# leaf_litte is % of leaf covering on the ground
# ground_veg is the undercover vegetation 1m down
# small_trun small large trees ...
# large_trun are the large trees in 2m radio from the frog
# fallen_tre are the numbre of fallen trees in 2m radio around the frog
# woody other type of dead wood in %
#Microhab_t the temp where the frog is standing
# Ambient_te the temp in the area
# canopy_cov % of canopy where the frog is standing
# Elevation general elevation of the location (gps)
# dist_to_ph distance to phytothelmata (water pools for babies)
# dist_to_tr distance to fallen tree ? 
# dust
# HFI human footprint index
# ind_per_km indivudals per km 


## subset and some data handling
PS_w <- subset_samples(PS, State=="Wild")
PS_w

PS_w@sam_data$Infected[61] <- "No"

PS_w@sam_data$Age <- "Adult"
PS_w@sam_data$Age[PS_w@sam_data$Sex=="Juvenile"] <- "Juvenile"
PS_w@sam_data$Age

PS_e <- subset_samples(PS, State=="Environment")
PS_e

PS_c <- subset_samples(PS, State=="Colony")
PS_c

## Data exploration and viz
# ordination plot
PS.ord <- vegan::metaMDS(PS@otu_table, distance="bray", try=250, k=3)
data.scores <- as.data.frame(vegan::scores(PS.ord)$sites)
df <- as.data.frame(PS@sam_data)

ordi_plot <- ggplot(data.scores, aes(NMDS1, NMDS2))+
    geom_point(size=2, alpha=0.8, aes(colour=df$State))+
    scale_color_manual(values=c("#002642", "#840032", "#e59500"))+
    theme_classic()+
    theme(legend.position="bottom")


# composition plot
taxW.df <- data.frame(tax_table(PS_w))
taxW.df %>% count(Phylum) -> taxW.df
taxW.df$Phylum[is.na(taxW.df$Phylum)] <- "Unknown_phylum"
taxW.df$Phylum <- as.factor(taxW.df$Phylum)

PS.t <- microbiome::transform(PS, "compositional")
ps.df <- psmelt(PS.t)

ps.df$Lab_ID <- as.factor(ps.df$Lab_ID)
ps.df$Phylum <- as.factor(ps.df$Phylum)
ps.df$State <- as.factor(ps.df$State)

all(levels(taxW.df$Phylum) == levels(ps.df$Phylum))

#ps.df$State[ps.df$Lab_ID2]

levels(ps.df$Phylum)# 54 phyla

(hcl.colors(54, palette = "Set 2"))

E1 <- ggplot(taxW.df, aes(x=Phylum, y=n, fill=Phylum))+
    geom_col()+
                                        #    scale_fill_paletteer_d("awtools::a_palette")+
    scale_fill_manual(values=hcl.colors(54, palette = "Set 2"))+
    coord_polar()+
    scale_y_log10("Number of cASVs per phylum")+
    theme_bw()+
    guides(fill="none")

E1

# compositin plot with relative abundances
ps.df <- tax_glom(PS, taxrank="Phylum")
ps.df <- transform(ps.df, "compositional") # transform to make all taxa from 0 to 1

ps.df <- psmelt(ps.df)
head(ps.df)
ps.df$Phylum <- as.factor(ps.df$Phylum)

E <- ggplot(ps.df, aes(x=Lab_ID, y=Abundance, fill=Phylum))+
    geom_bar(stat="identity")+
#    scale_fill_viridis(discrete=T)+
#    scale_fill_paletteer_d("ggsci::default_igv")+
    scale_fill_manual(values=hcl.colors(54, palette = "Set 2"))+
    labs(x="Samples", y="Relative abundance")+
    facet_grid(.~State, scales = "free_x", space="free_x")+
#    scale_x_discrete(labels = ps.df$State)+
    theme_classic()+
    guides(fill="none")+
    theme(axis.text.x = element_blank())

E

Fig1_1 <- plot_grid(ordi_plot, E1,labels=c("A", "B"))

Fig1 <- plot_grid(Fig1_1, E, ncol=1)

ggsave("fig/Figure1.pdf", Fig1, width=180, height=220, units="mm")



##### richness analyses"
########################
# first exploration between colony, wild and environment and coloured as infected and uninfected
p <- plot_richness(PS, x="State", measures=c("Chao1", "Shannon"))
p$layers <- p$layers[-1]
p+
    geom_boxplot(outlier.shape=NA)+
    geom_point(size=3, alpha=0.7, aes(color=Infected))+
    scale_fill_brewer(palette="Set1")+
    theme_classic() # no obvious differences, ignoring this for now and focusing only on the wild samples

# richness analysis for wild samples
rich.df <- estimate_richness(PS_w, split=TRUE, measures=c("Chao1", "Shannon"))
rich.df$HFI <- PS_w@sam_data$HFI
rich.df$Infect <- as.factor(PS_w@sam_data$Infected)
rich.df$InfL <- PS_w@sam_data$Bd.load_co
rich.df$Age <- as.factor(PS_w@sam_data$Age)
rich.df$BMI <- PS_w@sam_data$SMI.rob
rich.df$Loc <- as.factor(PS_w@sam_data$Site)

### the richness models

Chao1.rich <- lm(Chao1~Infect*HFI+Age+BMI+Loc, data=rich.df)
summary(Chao1.rich)
anova(Chao1.rich, test="LRT") # log likelihood ratio test 
plot(Chao1.rich) # looks good

Shannon.rich <- lm(Shannon~Infect*HFI+Age+BMI+Loc, data=rich.df)
summary(Shannon.rich)
anova(Shannon.rich, test="LRT")
plot(Shannon.rich) # also looks good


#########################################################################
##### Beta-diversity analysis
##### this is the distance based models. First for all, then only for adults
doData <- FALSE # leave as false to jump the heavy modeling

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

## 5) Age
AgeM<-array(0,c(nrow(metadt),nrow(metadt)))
for(i in 1:nrow(metadt)){
    for(j in 1:nrow(metadt)){
        if(metadt$Age[i]==metadt$Age[j]){
            AgeM[i,j]= "1"
        } else{
            AgeM[i,j]= "0"
        }
    }
}
AgeM<-c(AgeM[lower.tri(AgeM)])

## 6) HFI
hfiM <- c(dist(metadt$HFI))

##7) Parasite abundance
infM <- c(dist(metadt$Bd.load_co))

#Combine these vectors into a data frame
data.dyad<-data.frame(Jac=jac, Bray=bra, spatial=spa, BMI=bdm, Age=AgeM, HFI=hfiM, Infec=infM)

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

doModel <- FALSE
if(doModel){
# Now we model Presence/absence
modelJ<-brm(Jac~ 1+ spatial+ BMI + Age+ HFI+ Infec+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)

saveRDS(modelJ, "tmp/modelJ.rds")

# Now we model Abundances
modelB<-brm(Bray~ 1+ spatial+ BMI +Age+ HFI * Infec+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)

saveRDS(modelB, "tmp/modelB.rds")
} else{modelJ <- readRDS("tmp/modelJ.rds")
    modelB <- readRDS("tmp/modelB.rds")
}

###### models only for adults
doSex <- FALSE
PS_A <- subset_samples(PS_w, Age=="Adult")
summary(as.factor(PS_A@sam_data$Sex))
if(doSex){
metadt <- PS_A@sam_data
key <- sample_names(PS_A)

#1) Jaccard distances
JACM <- as.matrix(phyloseq::distance(PS_A, method="jaccard", type="samples", binary=T))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
# sanity check
all(rownames(JACM)==key)
jac <- c(as.dist(JACM))

#2) Bray distances
BRAM <- as.matrix(phyloseq::distance(PS_A, method="bray", type="samples"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
BRAM <- 1-BRAM
# sanity check
all(rownames(BRAM)==key)
bra <- c(as.dist(BRAM))

## 3) Spatial distance matrix
spa <- c(dist(metadt$Lon, metadt$Lat, method="euclidean"))

## 4) Bodymass distance matrix
bdm <- c(dist(metadt$SMI.rob))

## 5) Sex
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

# Now we model Presence/absence
modelJ_A<-brm(Jac~ 1+ spatial+ BMI + Sex+ HFI+ Infec+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)

saveRDS(modelJ_A, "tmp/modelJ_adults.rds")

# Now we model Abundances
modelB_A<-brm(Bray~ 1+ spatial+ BMI + Sex+ HFI * Infec+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)

saveRDS(modelB_A, "tmp/modelB_adults.rds")
} else{modelJ_A <- readRDS("tmp/modelJ_adults.rds")
    modelB_A <- readRDS("tmp/modelB_adults.rds")
}

modelB # this is bray-curtis distance model for all
modelB_A # this is bray-curtis distance model for only adults, to see sex effects

# ploting
abunM<-summary(modelB)$fixed
res.df<-data.frame(Effect=rownames(abunM),
                   Estimate=abunM$Estimate,
                   lCI=abunM$'l-95% CI',
                   uCI=abunM$'u-95% CI')
write.csv(res.df, "tmp/modelB.csv", row.names=FALSE) # saving model results

abunMA<-summary(modelB_A)$fixed
res.dfA<-data.frame(Effect=rownames(abunMA),
                   Estimate=abunMA$Estimate,
                   lCI=abunMA$'l-95% CI',
                   uCI=abunMA$'u-95% CI')
write.csv(res.dfA, "tmp/modelB_A.csv", row.names=FALSE) # saving model results

res.df <- res.df[-1,]


coul <- c("#591d44", "#9e217c", "#f28705", "#b5c516", "#f21f0c", "#e91e63")

Bray_eff <- ggplot(res.df, aes(x=Estimate, y=Effect, colour=Effect))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_point(size=3)+
    geom_errorbar(aes(xmin=lCI, xmax=uCI),
                  linewidth=1, width=0.4)+
             scale_colour_manual(values=coul)+
                labs(x="Estimate", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")
    

Bray_eff

ggplot2::ggsave(file="fig/Bray_eff.pdf", Bray_eff,  width = 190, height = 150, dpi = 300, units="mm")

### alternative vizualisation
modelB_transformed <- ggs(modelB)

cat <- filter(modelB_transformed, Parameter %in% c("b_spatial", "b_BMI", "b_Age1", "b_HFI", "b_Infec", "b_HFI:Infec"))

                                        # let's not forget to subset here
cat2 <- cat[cat$Iteration>1000,]

cat2$Parameter <- droplevels(cat2$Parameter)

unique(cat2$Parameter)


Fig2 <- ggplot(cat2, aes(x = value, y=Parameter, fill=Parameter))+
    geom_density_ridges(rel_min_height = 0.005, scale=5, alpha=0.8)+
#    facet_wrap(~Parameter)+
    geom_vline(xintercept = 0, col  = "black", size = 1, linetype="dashed")+
    geom_boxplot(outlier.shape = NA, width=0.2)+
    scale_fill_manual(values=coul)+
    xlab("Posterior probability distribution")+
    ylab("")+
    theme_classic()+
    theme(legend.position="none")

Fig2

###########################################################
## random forest regression for infection
###########################################################
PS.t <- transform(PS_w, "compositional") # transform to make all taxa from 0 to 1
otu.t <- PS.t@otu_table

colnames(otu.t) <- paste("ASV", seq(1, length(colnames(otu.t))), PS.t@tax_table[,6], sep="_")

df <- data.frame(otu.t,
                 InfC=as.factor(PS.t@sam_data$Infected),
                 BMI=PS.t@sam_data$SMI.rob,
                 Age=as.factor(PS.t@sam_data$Age),
                 HFI=PS.t@sam_data$HFI,
                 Coloration=as.factor(PS.t@sam_data$Coloration),
                 leaf_litte=PS.t@sam_data$leaf_litte,
                 ground_veg=PS.t@sam_data$ground_veg,
                 small_trun=PS.t@sam_data$small_trun,
                 large_trun=PS.t@sam_data$large_trun,
                 fallen_tre=PS.t@sam_data$fallen_tre,
                 woody=PS.t@sam_data$woody,
                 Microhab_t=PS.t@sam_data$Microhab_t,
                 canopy_cov=PS.t@sam_data$canopy_cov,
                 Elevation=PS.t@sam_data$Elevation,
                 dist_to_ph=PS.t@sam_data$dist_to_ph,
                 dist_to_tr=PS.t@sam_data$dist_to_tr,
                 ind_per_km=PS.t@sam_data$ind_per_km)

set.seed(123)
trainIndex <- createDataPartition(df$InfC, p=0.8, list=FALSE, times=1)
Inf_df_train <- df[trainIndex,]
Inf_df_test <- df[-trainIndex,]

tgrid <- expand.grid(
    mtry = 5:100)

set.seed(123)
fitControl <- trainControl(method="cv", number=5, search="grid")

doRF <- FALSE

fitControl <- trainControl(method="cv", number=5, classProbs=TRUE,
                           summaryFunction=twoClassSummary,
                           savePredictions="final",
                           search="grid")

#devtools::install_github('topepo/caret/pkg/caret')

if (doRF){
    set.seed(1111)
    rfFit2 <- train(InfC~., data=Inf_df_train,
                    method="ranger",
                    metric="ROC",
                    trControl=fitControl,
#                    tuneGrid=tgrid,
                    tuneGrid=expand.grid(
                        mtry=c(5, 10, 20, 30, 40, 50, 100),
                        splitrule="gini",
                        min.node.size=c(1,5, 10, 20)),
                    importance="permutation")
    saveRDS(rfFit2, "tmp/rfFit2.rds")
} else
rfFit2 <- readRDS("tmp/rfFit2.rds")


print(rfFit2)

#print(rfFit1$finalModel)

test_predictions <- predict(rfFit2, newdata=Inf_df_test)

print(postResample(test_predictions, Inf_df_test$InfC))

confusionMatrix(test_predictions, Inf_df_test$InfC)


### let's see the 20 featured that are most associated with infection
imp <- varImp(rfFit1)
imp_df <- imp$importance
imp_df$feature <- rownames(imp_df)
rownames(imp_df) <- NULL
all(imp_df$No==imp_df$Yes) # sanity check

imp_df$No <- NULL
imp_df <- imp_df[order(-imp_df$Yes),]
imp20 <- imp_df[1:20,]
imp20 <- droplevels(imp20)

imp20

#plot(varImp(rfFit2), top = 20)

Taxa20 <- data.frame(seq=taxa_names(PS.t)[as.numeric(rownames(imp20))], name=imp20$feature, importance=imp20$Yes)

library(seqinr)
write.csv2(Taxa20, "tmp/Taxatop20.csv")
write.fasta(as.list(Taxa20$seq), Taxa20$name, "tmp/Taxatop20.fasta")

topImp <-
    ggplot(imp20, aes(y=reorder(feature, Yes), x=Yes))+
    geom_segment( aes(yend=feature, xend=0)) +
    geom_point(size=4, color="orange")+
    labs(x="feature importance", y="")+
    theme_classic()

topImp

