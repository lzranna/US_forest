library(phyloseq)
library(bendiR)
library(ggplot2)
library(ggpubr)
library(readxl)
library(vegan)
library(viridis)
library(dplyr)
library(performance)
library(RColorBrewer)
library(phylogeo)
library(ggordiplots)
library(lme4)
library(sjPlot)
library(microbiome)
library(eulerr)
library(phylosmith)
library(tidyverse)
library(microViz)
library(MuMIn)
library(jtools)

setwd("/Users/u2061312/Documents/FOREST_Ryan/R")
load("/Users/u2061312/Documents/FOREST_Ryan/R/FOREST_phyloseq.R")
ps
ps_RA <- transform_sample_counts(ps, function(x) x/sum(x)*100)
get_taxa_unique(ps, "Genus")


ECM_ra <- subset_taxa(ps_RA, Genus=="g__Cenococcum" | Family=="f__Elaphomycetaceae" |
                     Genus=="g__Hydnotrya" | Genus=="g__Hygrophorus" |
                     Genus=="g__Helvella" | Genus=="g__Hebeloma" |
                     Genus=="g__Delastria" |Genus=="g__Hymenogaster" |
                     Genus=="g__Hydnobolites" | Genus=="g__Phaeocollybia" |
                     Genus=="g__Scabropezia" | Family=="f__Inocybaceae" |
                     Genus=="g__Genea" | Genus=="g__Lyophyllum" |
                     Genus=="g__Otidea" | Genus=="g__Stephanospora" |
                     Genus=="g__Pulvinula" | Genus=="g__Piloderma" |
                     Genus=="g__Sphaerosporella" | Genus=="g__Amphinema" |
                     Genus=="g__Tarzetta" | Genus=="g__Tylospora" |
                     Genus=="g__Tuber" | Genus=="g__Astraeus" |
                     Genus=="g__Phlegmacium" | Family=="f__Boletaceae" |
                     Genus=="g__Laccaria" | Genus=="g__Gyroporus" |
                     Genus=="g__Melanogaster" | Genus=="g__Paxillus" |
                     Genus=="g__Cantharellus" | Genus=="g__Scleroderma" |
                     Family=="f__Rhizopogof" | Genus=="g__Hydnum" |
                     Genus=="g__Suillus" | Genus=="g__Epulorhiza" |
                     Genus=="g__Sclerogaster" | Family=="f__Gautieriaceae" |
                     Genus=="g__Coltricia" | Genus=="g__Hysterangium" |
                     Genus=="g__Leucogaster" | Family=="f__Russulaceae" |
                     Genus=="g__Tremellodendron" | Family=="f__Bankeraceae" |
                     Family=="f__Thelephoraceae" | Genus=="g__Endogone")
ECM <- subset_taxa(ps, Genus=="g__Cenococcum" | Family=="f__Elaphomycetaceae" |
                        Genus=="g__Hydnotrya" | Genus=="g__Hygrophorus" |
                        Genus=="g__Helvella" | Genus=="g__Hebeloma" |
                        Genus=="g__Delastria" |Genus=="g__Hymenogaster" |
                        Genus=="g__Hydnobolites" | Genus=="g__Phaeocollybia" |
                        Genus=="g__Scabropezia" | Family=="f__Inocybaceae" |
                        Genus=="g__Genea" | Genus=="g__Lyophyllum" |
                        Genus=="g__Otidea" | Genus=="g__Stephanospora" |
                        Genus=="g__Pulvinula" | Genus=="g__Piloderma" |
                        Genus=="g__Sphaerosporella" | Genus=="g__Amphinema" |
                        Genus=="g__Tarzetta" | Genus=="g__Tylospora" |
                        Genus=="g__Tuber" | Genus=="g__Astraeus" |
                        Genus=="g__Phlegmacium" | Family=="f__Boletaceae" |
                        Genus=="g__Laccaria" | Genus=="g__Gyroporus" |
                        Genus=="g__Melanogaster" | Genus=="g__Paxillus" |
                        Genus=="g__Cantharellus" | Genus=="g__Scleroderma" |
                        Family=="f__Rhizopogof" | Genus=="g__Hydnum" |
                        Genus=="g__Suillus" | Genus=="g__Epulorhiza" |
                        Genus=="g__Sclerogaster" | Family=="f__Gautieriaceae" |
                        Genus=="g__Coltricia" | Genus=="g__Hysterangium" |
                        Genus=="g__Leucogaster" | Family=="f__Russulaceae" |
                        Genus=="g__Tremellodendron" | Family=="f__Bankeraceae" |
                        Family=="f__Thelephoraceae" | Genus=="g__Endogone")                   


AM_ra <- subset_taxa(ps_RA, Phylum=="p__Glomeromycota")
AM <- subset_taxa(ps, Phylum=="p__Glomeromycota")


#this needs fixing it has the wrong OTU table and tax table see bellow
load("/Users/u2061312/Documents/FOREST_Ryan/R/group_Tetracladium_ps.R")
group 

load("/Users/u2061312/Documents/FOREST_Ryan/R/group_Tetracladium_ps_RA.R")
group_RA
grop_tt <- tax_table(group_RA)
group_meta <- sample_data(group)
OTUdf2

#rename tetracladium asvs
name_Tetra <- Tetracladium_ps
taxa_names(name_Tetra) <- paste0("ITS_ASV_MANE_", seq(ntaxa(name_Tetra)))
otu_table(name_Tetra, taxa_are_rows = F)

#rename all asvs for deposit
otu_table(ps, taxa_are_rows = F)
fungi_rest = subset_taxa(ps, Genus!="g__Tetracladium"| is.na(Genus))
tetra_rest = subset_taxa(ps, Genus=="g__Tetracladium")
fungi_rest_name<- fungi_rest
taxa_names(fungi_rest_name) <- paste0("ITS_ASV_MANE_FUNGI_", seq(ntaxa(fungi_rest_name)))
otu_table(fungi_rest_name, taxa_are_rows = F)

phylo_deposit <- merge_phyloseq(name_Tetra, fungi_rest_name)
phylo_fasta <- merge_phyloseq(tetra_rest, fungi_rest)

OTU_depo = as(otu_table(phylo_fasta), "matrix")
if(taxa_are_rows(phylo_fasta)){OTU_depo <- t(OTU_depo)}
OTUdf1 = as.data.frame(OTU_depo)
Sample1 <- rownames(OTUdf1)
rownames(OTUdf1) <- NULL
OTUdf21 <- cbind(Sample1,OTUdf1)

write.csv(OTUdf21, "OTUtable_fasta.csv")

#new group
otu_table <- otu_table(OTU1, taxa_are_rows = F)

new_group <- phyloseq(otu_table, group_meta, grop_tt)
  
sample_variables(ps)
get_taxa_unique(ps, "Phylum")

Tetracladium_ps <- subset_taxa(ps, Genus=="g__Tetracladium")

Tetracladium_ps_RA <- subset_taxa(ps_RA, Genus=="g__Tetracladium")

UpMin <- subset_samples(Tetracladium_ps_RA, depth ==5)
metadata_um <- as(sample_data(UpMin), "data.frame")
UpMin_n <- subset_samples(Tetracladium_ps, depth ==5)


LowMin <- subset_samples(Tetracladium_ps_RA, depth ==15)
metadata_lm <- as(sample_data(LowMin), "data.frame")

Enzime <- subset_samples(Tetracladium_ps_RA, BG_Fall !="NA")
metadata_e <- as(sample_data(Enzime), "data.frame")
Enzime_D <- subset_samples(Tetracladium_ps, BG_Fall !="NA")


#Remove samples with no metadata
sample_variables(Tetracladium_ps_RA)
phylo <- subset_samples(Tetracladium_ps_RA, horizon !="")
phylo <- subset_samples(phylo, myco !="")
phylo <- subset_samples(phylo, ecm_tree.percent !="")
phylo <- subset_samples(phylo, n.percent !="")
phylo <- subset_samples(phylo, c.percent !="")
phylo <- subset_samples(phylo, s.percent !="")
phylo <- subset_samples(phylo, cn !="")
phylo <- subset_samples(phylo, cs !="")
phylo <- subset_samples(phylo, ns !="")
phylo <- subset_samples(phylo, d15N !="")
phylo <- subset_samples(phylo, d13C !="")
phylo <- subset_samples(phylo, d34S !="")
phylo <- subset_samples(phylo, ph !="")

#REMOVE taxa and samples with 0 no reads  
phylo <- prune_taxa(taxa_sums(phylo) > 0, phylo)
phylo <- prune_samples(sample_sums(phylo)>0,phylo)


#Make sure this phyloseq matches the phyloseq you are using.
metadata <- as(sample_data(phylo), "data.frame")

#rarefy
## Rarefaction curve, ggplot style
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed.
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}


plot_rare <- ggrare(ps, step = 100, color = "site", label = NULL, se = FALSE)

ggsave(filename = "RareCurve.pdf", plot = plot_rare,
       scale = 1,
       dpi = 300)


#write OTU ra table
OTU1 = as(otu_table(group_RA), "matrix")
if(taxa_are_rows(group_RA)){OTU1 <- t(OTU1)}
OTUdf = as.data.frame(OTU1)
Sample <- rownames(OTUdf)
rownames(OTUdf) <- NULL
OTUdf2 <- cbind(Sample,OTUdf)

write.csv(OTUdf2, "OTUtable_RA.csv")

###################
#Richness
diversity_plot4<-function(physeq,grouping,diversity=c("fishers","shannon","simpson","observed"),test=c("anova","kruskal"),correction_method=stats::p.adjust.methods){
  if(!phyloseq::taxa_are_rows(physeq)){
    phyloseq::otu_table(physeq)<-t(phyloseq::otu_table(physeq))
  }
  test <- base::match.arg(test,c("anova","kruskal"))
  diversity <- base::match.arg(diversity,c("fishers","shannon","simpson","observed"))
  if(length(correction_method)>1){
    corr_specified<-FALSE
  } else{
    corr_specified<-TRUE
  }
  if(length(test)==2){
    warning("Defaulting to ANOVA")
    test<-"anova"
  }
  if(!length(grep(grouping,colnames(phyloseq::sample_data(physeq))))==0){
    if(diversity=="fishers"){
      div<-vegan::fisher.alpha(t(phyloseq::otu_table(physeq)))
    }  else if(diversity=="observed"){
      div<-colSums(otu_table(physeq) != 0)
    } else {
      div<-vegan::diversity(t(phyloseq::otu_table(physeq)),diversity)
    } 
    if(test=="anova"){
      anova<-stats::aov(div ~ as.character(phyloseq::sample_data(physeq)[[grouping]]))
      post_hoc<-stats::TukeyHSD(anova)
      means<-stats::aggregate(anova$model[, 1], list(anova$model[,2]), mean)
      letters<-compact_letters(phyloseq::sample_data(physeq)[[grouping]],post_hoc)
    } else{
      correction_method <- base::match.arg(correction_method)
      if(corr_specified == FALSE){
        warning("Multiple testing correction post-hoc p-values with FDR")
        correction_method<-"fdr"
      }
      kruskal<-stats::kruskal.test(div ~ as.factor(phyloseq::sample_data(physeq)[[grouping]]))
      post_hoc<-DescTools::DunnTest(div~as.factor(phyloseq::sample_data(physeq)[[grouping]]),method=correction_method)
      diversity_table<-cbind.data.frame(div,group=as.character(phyloseq::sample_data(physeq)[[grouping]]))
      means<-stats::aggregate(diversity_table[,"div"], list(diversity_table[,"group"]), mean)
      letters<-compact_letters(phyloseq::sample_data(physeq)[[grouping]],post_hoc)
    }
    
    names(means)<-c("Group","mean")
    
  } else{
    stop(paste("Grouping factor entered does not exist in phyloseq::sample_data. Check phyloseq::sample_data(",substitute(physeq),").",sep=""))
  }
  if(!exists("diversity_table")){
    diversity_table<-cbind.data.frame(div,group=as.character(phyloseq::sample_data(physeq)[[grouping]]))
  }
  
  if(typeof(letters)=="character"){
    letters<-cbind.data.frame(sample=names(letters),letters)
  } else{
    letters<-cbind.data.frame(sample=names(letters$Letters),letters=letters$Letters)
    print(letters)
  }
  plot<-ggplot2::ggplot(data=diversity_table)+ggplot2::geom_boxplot(aes(x=group,y=div,fill=group, colour=group),alpha=0.3,outlier.shape=NA)+ggplot2::geom_jitter(aes(x=group,y=div, colour=group),alpha=0.5)+ggplot2::geom_text(data=letters,aes(x=sample,y=Inf,label=letters,vjust = 1,fontface="bold"),size=4)+theme_journal()+scale_y_continuous("Observed species")+ theme(legend.title = element_blank())#+
  #theme(axis.title.x=element_blank(),
  #axis.text.x=element_blank(),
  #axis.ticks.x=element_blank())
  return(list(POSTHOC=post_hoc,PLOT=plot))
}

alpha <- diversity_plot4(UpMin,"site",diversity="observed",test="kruskal", correction_method = "none") 
plot1 <- alpha$PLOT
plot2 <- plot1  +
  theme_pubr() +
  theme(legend.position="right", axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(), legend.title = element_blank())

alpha$POSTHOC
ggsave(filename = "Richness_site.pdf", plot = plot2,
       scale = 1.5,
       dpi = 300)

beta <- diversity_plot4(Tetracladium_ps_RA,"myco",diversity="observed",test="kruskal", correction_method = "none") 
plot1b <- beta$PLOT
plot2b <- plot1b  +
  theme_pubr() +
  theme(legend.position="right", axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(), legend.title = element_blank())

beta$POSTHOC
ggsave(filename = "Richness_myco.pdf", plot = plot2b,
       scale = 1.5,
       dpi = 300)


##################
#NMDS
set.seed(12345)
GP.ord2 <- ordinate(phylo, distance="bray", method="MDS")
p2 = plot_ordination(phylo, GP.ord2, type="samples", color="site",label=NULL,axes=c(1, 2))+
  geom_point(aes(color=site),size=2) +
  ggtitle("PCoA")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_ellipse(geom = "polygon", type = "t", alpha = 0.3, aes(fill = site))+
  geom_point(size=5)

p2

library(pairwiseAdonis)
ps_tr <- microbiome::transform(phylo, "clr")

ps_dist_matrix <- phyloseq::distance(phylo, method ="bray")

vegan::adonis2(ps_dist_matrix ~ phyloseq::sample_data(ps_tr)$site)

site_pw <- pairwise.adonis(ps_dist_matrix, phyloseq::sample_data(ps_tr)$site)
write_csv(site_pw, "site_pw.csv")




GP.ord3 <- ordinate(UpMin_full, distance="bray", method="MDS")

p3 = plot_ordination(UpMin_full, GP.ord3, type="samples", color="SoilT",shape="site",label=NULL,axes=c(1, 2))+
  geom_point(aes(color=SoilT),size=5) +
  ggtitle("PCoA")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_shape_manual(values=c(3, 15, 16, 17,18))

p3

p4 = plot_ordination(phylo, GP.ord2, type="samples", color="myco",label=NULL,axes=c(1, 2))+
  geom_point(aes(color=myco),size=5) +
  ggtitle("PCoA")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_ellipse(geom = "polygon", type = "t", alpha = 0.3, aes(fill = myco))
  

p4


vegan::adonis2(ps_dist_matrix ~ phyloseq::sample_data(ps_tr)$myco)

myco_pw <- pairwise.adonis(ps_dist_matrix, phyloseq::sample_data(ps_tr)$myco)
write_csv(myco_pw, "myco_pw.csv")

ggsave(filename = "pcoa_site.pdf", plot = p2,
       scale = 1.5,
       dpi = 300)

ggsave(filename = "pcoa_myco.pdf", plot = p4,
       scale = 1.5,
       dpi = 300)

#PERMANOVA use cleaned RA dataset (NO NAs) this case phylo

#Run this function to convert physeq to vegan format
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

sample_variables(phylo)


#Test each variable (categorical or continuous) separately,

adonis2(veganotu(phylo) ~ depth, data = metadata, method="bray", permutations = 100)
adonis2(veganotu(phylo) ~ horizon, data = metadata, method="bray", permutations = 100)
adonis2(veganotu(phylo) ~ site, data = metadata, method="bray", permutations = 100)
adonis2(veganotu(phylo) ~ myco, data = metadata, method="bray", permutations = 100)
adonis2(veganotu(phylo) ~ ecm_tree.percent, data = metadata, method="bray", permutations = 100)
adonis2(veganotu(phylo) ~ n.percent, data = metadata, method="bray", permutations = 100)
adonis2(veganotu(phylo) ~ c.percent, data = metadata, method="bray", permutations = 100)
adonis2(veganotu(phylo) ~ s.percent, data = metadata, method="bray", permutations = 100)
adonis2(veganotu(phylo) ~ ph, data = metadata, method="bray", permutations = 100)


#Then put together in order of highest R first.

set.seed(12345)

adonis <- adonis2(veganotu(phylo) ~ myco+ecm_tree.percent+n.percent+ph+
                    c.percent+s.percent+depth+horizon,
                 data = metadata, method="bray", permutations = 100)


adonis

#Adjusted p values. Not sure if needed. n = the number of variables added. Can change to holm, fdr etc. Canoco uses FDR. 
adonis$aov.tab$`Adj_P` <- p.adjust(adonis$aov.tab$`Pr(>F)`, method = 'fdr', n = 9)
#See table with adjusted p values for multiple comparisons included. 

adonis

#Write output to text.

out <- capture.output(adonis)
cat(out, file="adonis_no_site.txt", sep="\n")

#adonis with upmin (0-5) only.


sample_variables(UpMin)
UpMin_full <- subset_samples(UpMin, horizon !="")
UpMin_full <- subset_samples(UpMin_full, myco !="")
UpMin_full <- subset_samples(UpMin_full, ecm_tree.percent !="")
UpMin_full <- subset_samples(UpMin_full, n.percent !="")
UpMin_full <- subset_samples(UpMin_full, c.percent !="")
UpMin_full <- subset_samples(UpMin_full, s.percent !="")
UpMin_full <- subset_samples(UpMin_full, cn !="")
UpMin_full <- subset_samples(UpMin_full, cs !="")
UpMin_full <- subset_samples(UpMin_full, ns !="")
UpMin_full <- subset_samples(UpMin_full, d15N !="")
UpMin_full <- subset_samples(UpMin_full, d13C !="")
UpMin_full <- subset_samples(UpMin_full, d34S !="")
UpMin_full <- subset_samples(UpMin_full, ph !="")

#REMOVE taxa and samples with 0 no reads  
UpMin_full <- prune_taxa(taxa_sums(UpMin_full) > 0, UpMin_full)
UpMin_full <- prune_samples(sample_sums(UpMin_full)>0,UpMin_full)


#Make sure this phyloseq matches the phyloseq you are using.
metadata_upmin <- as(sample_data(UpMin_full), "data.frame")
metadata_upmin <- select(metadata_upmin, -c(NAG_Spring,NAG_Summer,NAG_Fall,BG_Summer,BG_Spring,
                                            BG_Fall,Phenox_Spring,Phenox_Summer,Phenox_Fall,
                                            Pheno_Spring,Pheno_Summer,Pheno_Fall, RootN,
                                            RootC,RootCN, AP_Spring, AP_Summer, AP_Fall))



adonis2(veganotu(UpMin_full) ~ site, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ myco, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ ecm_tree.percent, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ n.percent, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ c.percent, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ s.percent, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ ph, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ LitCN, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ LitLCI, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ LitLigN, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Sand, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Silt, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Clay, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Feo, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Alo, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Gwc, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Whc, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Nmin, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Nit, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ SoilT, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Glu, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Gal, data = metadata_upmin, method="bray", permutations = 100)
adonis2(veganotu(UpMin_full) ~ Mur, data = metadata_upmin, method="bray", permutations = 100)




adonis_upmin <- adonis2(veganotu(UpMin_full) ~site+SoilT+Whc+Glu+Gal+ myco+ Mur+ Clay+ Gwc+LitLCI+Sand+Nmin+Alo+
                    n.percent+ph+ecm_tree.percent+LitLigN+s.percent,
                  data = metadata_upmin, method="bray", permutations = 100)


out <- capture.output(adonis_upmin)
cat(out, file="adonis_upmin.txt", sep="\n")


#################
#RDA analyses source: http://www.hiercourse.com/docs/microbial/04_betaDiversity_multTables.html
# extract the OTU table
Y <- veganotu(UpMin)

# extract the sample data from the 'phyloseq' object
# then remove the 'SampleID' column
x <- select(mod_data_norm, -c(x))

# RDA
res <- rda(Y ~ ., data=x)
resplot1 <- plot(res, xlab=ord_labels(res)[1], ylab=ord_labels(res)[2])
summary(res, display=NULL)

#There are multiple approaches to select only the variables that are explaining variation efficiently.
# 1.calculate variance inflation
# variables with scores >10 are redundant
sort(vif.cca(res))

# 2. calculate fit - Fitting environmental vectors/factors onto an ordination 
#envfit(Y ~ ., data=X)

# 3. stepwise selection - set up full and null models for 'ordistep'
# full model
rda1 <- rda(Y ~ ., data=x)

# intercept-only (null) model
rda0 <- rda(Y ~ 1, data=x)

# perform forward and backward selection of explanatory variables
# output not shown
step.env <- ordistep(rda0, scope=formula(rda1), direction='both')
# look at the significant variables 
step.env$anova
# code to get variable names from 'ordistep' and 'envfit' results
vars.ordistep <- gsub('^. ', '', rownames(step.env$anova))

#vars.envfit <- names(which(vif.cca(res) <= 10))

#vars <- unique(c(vars.ordistep, vars.envfit))

# select variables to keep from table 'Y'
#X1 <- X[, vars]
#str(X1)
#RDA REDUCED VARIABLES
#res <- rda(Y ~ ., data=X1)

# summary of the results
#summary(res, display=NULL)
#anova(res)


# set up dataframes for plotting the results
#sit <- cbind(dat, scores(res, display='sites'))
#spp <- cbind(data.frame(tax_table(phylo)), scores(res, display='species'))
#vec <- data.frame(scores(res, display='bp'))

# use these to adjust length of arrows and position of arrow labels
#adj.vec <- 2
#adj.txt <- 2.5

# 'site' ordination
# p1 <- ggplot(sit, aes(x=RDA1, y=RDA2, color=c.percent, shape=site)) +
#   geom_point(size=3) +
#   scale_shape_manual(values=c(15,4,17,18,19,2,8,11))+
#   geom_segment(data=vec, inherit.aes=F, 
#                mapping=aes(x=0, y=0, xend=adj.vec*RDA1, yend=adj.vec*RDA2), 
#                arrow=arrow(length=unit(0.2, 'cm'))) + 
#   geom_text(data=vec, inherit.aes=F, 
#             mapping=aes(x=adj.txt*RDA1, y=adj.txt*RDA2, 
#                         label=c('ecm_tree.percent', 'c.percent ',"s.percent", 'ph'))) +
#   theme_bw()
# p1$labels$x <- ord_labels(res)[1]
# p1$labels$y <- ord_labels(res)[2]

#ggsave("finalRDA_ordi.pdf", p1, dpi = 300)

#Variation partitioning - more than two tables

# define the partitions
# soil variables from ordistep and envfit analyses
# X1 <- X[, vars]
# X1 <- select(X1, -ecm_tree.percent)
# habitat variables
#X2 <- select(dat, site)
# myco variables
#X3 <- select(dat, ecm_tree.percent)
# profile variables
#X4 <- select(dat, horizon)

# partition variation and plot the result
#vp <- varpart(Y, X1, X2, X3, X4)
#plot(vp, Xnames=c('Soil prop', 'Site', 'ECM percent', "Horizon"))

#We can test the significance of each partition by performing RDA and conditioning the OTU matrix with all of the other partitions than the one that we are focussing on in that analysis. For example:
# Test the 'Habitat' partition (X2)
# first argument contains the community data 
# second argument is the partition of interest (constraint)
# third argument is a column-bound dataframe of remaining partitions (condition)
#rda.x2 <- rda(Y, X2, cbind(X1, X3))
#rda.x1 <- rda(Y, X1, cbind(X2, X3))
#rda.x3 <- rda(Y, X3, cbind(X2, X1))

# test significance
# anova(rda.x2)
# anova(rda.x1)
# anova(rda.x3)
################
###########################MODELING
#correlation matrix RA upmin
Y_cor <- veganotu(UpMin)
Ysum_cor <- rowSums(Y_cor,na.rm = FALSE)

dat_cor <- data.frame(sample_data(UpMin))
dat_cor <- select(dat_cor, -cn)
dat_cor <- select(dat_cor, -cs)
dat_cor <- select(dat_cor, -ns)
dat_cor <- select(dat_cor, -d13C )
dat_cor <- select(dat_cor, -d15N )
dat_cor <- select(dat_cor, -d34S)
dat_cor <- select(dat_cor, -depth)

mod_data_cor <- merge(Ysum_cor, dat_cor, by ="row.names", all = TRUE )
rownames(mod_data_cor) <- mod_data_cor[,1]
mod_data_cor[,1] <- NULL
mod_data_cor

#grab data (has to be data frame)
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#select variables to keep (no NAs)
mod_data_cor2 = subset(mod_data_cor, select = c(x,ecm_tree.percent, n.percent, c.percent,
                                                s.percent,ph,LitCN,LitLCI,LitLigN,Sand,Silt,Clay,
                                                Feo,Alo, Gwc,Whc,Nmin,Nit,SoilT,Glu,Gal,Mur))


mod_data_norm <- as.data.frame(lapply(mod_data_cor2[1:22], min_max_norm))

res2 <- rcorr(as.matrix(mod_data_norm))

cp1<- ggcorrplot(res2$r, hc.order = TRUE, type = "lower",
                 p.mat = res2$P,
                 colors = c("#327C87", "white", "#992B19"))
cp1


ggsave(plot =cp1, file = "corr_ra.pdf" )

#add back site column
mod_data_norm2<-mod_data_norm
mod_data_norm2$site <- mod_data_cor$site


#site as random variable
mod1 <- lmer(x ~ecm_tree.percent+ n.percent+ c.percent+
             s.percent+ph+LitCN+LitLCI+LitLigN+Sand+Silt+
             Feo+Alo+ Gwc+Whc+Nmin+Nit+SoilT+Glu+Gal+Mur + (1|site), data = mod_data_norm2)

mod_data_norm2$CN <- mod_data_norm2$c.percent/mod_data_norm2$n.percent

mod_data_norm2 <- mod_data_norm2 %>% 
  filter_all(all_vars(!is.infinite(.)))

mod2 <- lm(x ~ecm_tree.percent+ CN+
               s.percent+ph+LitCN+LitLCI+LitLigN+Sand+Silt+
               Feo+Alo+ Gwc+Whc+Nmin+Nit+SoilT+Glu+Gal+Mur , data = mod_data_norm2)

mod3 <- lmer(x ~SoilT +(1|site), data = mod_data_norm2)

tab_model(mod1, show.est=TRUE)
plot_mod1 <- plot_model(mod1, sort.est=FALSE, type = "std", show.values=TRUE, value.size=2, value.offset=0.4)
plot_mod1

summary(mod2)

tab_model(mod3, show.est=TRUE)
plot_mod3 <- plot_model(mod3, sort.est=FALSE, type = "std", show.values=TRUE, value.size=2, value.offset=0.4)
plot_mod3


compare_performance(mod1, mod2, mod3, verbose = FALSE, rank = TRUE)


# Model_pSEM <- psem(
#   lm(x ~ ecm_tree.percent + c.percent + Sand + ph + LitLCI, data =mod_data_norm2),
#   lm(ph ~ ecm_tree.percent +c.percent +Sand+LitLCI, data =mod_data_norm2))
# 
# PSEM <- plot(
#   Model_pSEM,
#   return = F,
#   node_attrs = data.frame(shape = "rectangle", style = "rounded", fixedsize = F, color = "#2F4F4F", fillcolor = "white", fontsize = 8, penwidth = 2),
#   edge_attrs = data.frame(style = "solid", color = "black", arrowhead = "none"),
#   ns_dashed = T,
#   alpha = 0.05,
#   show = "std",
#   digits = 3,
#   add_edge_label_spaces = T,
#   title = "RA PSEM")
# 
# PSEM
# summary(Model_pSEM, test.statistic = "T")



###############
#model with richness
###########################MODELING
#grab data (has to be data frame)
observed_D_samples <- estimate_richness(UpMin_n, split = TRUE, measures = "Observed")
#observed_D_samples <- filter(observed_D_samples, Observed > 0)
tryaa<-cbind(observed_D_samples,dat_cor)

#select variables to keep (no NAs)
tryaa2 = subset(tryaa, select = c(Observed,ecm_tree.percent, n.percent, c.percent,
                                                s.percent,ph,LitCN,LitLCI,LitLigN,Sand,Silt,Clay,
                                                Feo,Alo, Gwc,Whc,Nmin,Nit,SoilT,Glu,Gal,Mur))


tryaa_norm <- as.data.frame(lapply(tryaa2[1:22], min_max_norm))

#add back site column
tryaa_norm2<-tryaa_norm
tryaa_norm2$site <- mod_data_cor$site

#Site as random variable

mod3R <- lmer(Observed ~  ecm_tree.percent+ n.percent+ c.percent+
                s.percent+ph+LitCN+LitLCI+LitLigN+Sand+Silt+
                Feo+Alo+ Gwc+Whc+Nmin+Nit+SoilT+Glu+Gal+Mur  + (1|site), data = tryaa_norm2)
tab_model(mod3R, show.est=TRUE)

plot_mod3R <- plot_model(mod3R, sort.est=FALSE, type = "std", show.values=TRUE, value.size=2, value.offset=0.4)
plot_mod3R

tryaa_norm2$CN <- tryaa_norm2$c.percent/tryaa_norm2$n.percent

tryaa_norm2 <- tryaa_norm2 %>% 
  filter_all(all_vars(!is.infinite(.)))


mod2R <- lm(Observed ~ecm_tree.percent+ CN +
             s.percent+ph+LitCN+LitLCI+LitLigN+Sand+Silt+
             Feo+Alo+ Gwc+Whc+Nmin+Nit+SoilT+Glu+Gal+Mur , data = tryaa_norm2)

summary(mod2R)

mod4R <- lm(Observed ~ph+LitLCI+LitLigN+Sand+
              Alo+SoilT , data = tryaa_norm2)

summary(mod4R)


compare_performance(mod3R, mod2R,mod4R, verbose = FALSE, rank = TRUE)
#stepwise variable selection
#  prevent fitting sub-models to different datasets
options(na.action = "na.fail")
dd<-dredge(mod2, trace = 2, extra = c("R^2", F = function(x)
  summary(x)$fstatistic[[1]]))

subset(dd, delta < 4)

# Visualize the model selection table:
par(mar = c(3,5,6,4))
plot(dd, labAsExpr = TRUE)
# Model average models with delta AICc < 4
model.avg(dd, subset = delta < 4)

#or as a 95% confidence set:
model.avg(dd, subset = cumsum(weight) <= .95) # get averaged coefficients

#'Best' model
summary(get.models(dd, 1)[[1]])


best_mod_ra <- lm(formula = x ~ Gal + LitCN + ph + s.percent + Silt + SoilT + 
                    1, data = mod_data_norm2)

summary(best_mod_ra)

dd_r<-dredge(mod2R, trace = 2, extra = c("R^2", F = function(x)
  summary(x)$fstatistic[[1]]))
summary(get.models(dd_r, 1)[[1]])

best_mod_obs <- lm(formula = Observed ~  Glu + LitCN + LitLCI + LitLigN + 
                     Nit + ph + Silt + Whc + 1, data = tryaa_norm2)
summary(best_mod_obs)

CMplotR <- check_model(
  best_mod_obs,
  dot_size = 2,
  line_size = 0.8,
  panel = TRUE,
  check = "all",
  alpha = 0.2,
  dot_alpha = 0.8,
  colors = c("#3aaf85", "#1b6ca8", "#cd201f"),
  theme = "see::theme_lucid",
  detrend = FALSE,
  verbose = TRUE)

plot_summs(best_mod_obs)
plot_summs(best_mod_ra)
all_lin <- plot_summs(best_mod_ra, best_mod_obs,
           scale = TRUE, robust = TRUE,
           model.names = c("Observed","RA"))

ggsave("both_lin_sum.pdf", all_lin, dpi = 300)



# model_coefsR <- coef(mod3R)$ site %>% 
#   rename(Intercept = `(Intercept)`, Slope = ecm_tree.percent) %>% 
#   rownames_to_column("site")
# 
# sleep_groups_raniR <- left_join(tryaa, model_coefsR, by = "site")
# 
# model_coef_plotR <- ggplot(data = sleep_groups_raniR, 
#                            mapping = aes(x = ecm_tree.percent, 
#                                          y = Observed, 
#                                          colour = site)) +
#   geom_point(na.rm = T, alpha = 0.5) +
#   geom_abline(aes(intercept = Intercept, 
#                   slope = Slope,
#                   colour = site),size = 1.5) +
#   theme_pubr()


ggsave("model_coef_plotR.pdf", model_coef_plotR, dpi = 300)
##############
#STACKED BAR PLOTS
phyla<- aggregate_taxa(group_RA, level = "Species")

p_reads <- plot_composition(phyla,
                            average_by = "site",
                            transform = "compositional") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))
p_reads

#Make barplot

#Merge by phyla, change Phylum to Class (or different ta level) if needed
phyla<-(tax_glom(new_group,"Species"))

#Merge by Sample type (column heading of your metadata). Don't worry about NAs.
merged_phyla<-merge_samples(phyla, "site")


otu_table(merged_phyla)<-otu_table((otu_table(merged_phyla))/(data.frame(sample_data(new_group)) %>% group_by(site) %>% tally())$n,taxa_are_rows=FALSE)

#In this one you have the choice of leaving all taxa that are above x% in at least one sample in, or just merging everything low abundance on a sample by sample basis, keep_all=TRUE

#merged_phyla<-merge_low_abun_taxa(merged_phyla,1,keep_all=TRUE)

#Melt for use with ggplot
melted_phyla<-psmelt(merged_phyla)

plot_phyla<-ggplot(data=melted_phyla, 
                   aes(fill=Species,x=Sample, y=Abundance))+
  geom_bar(stat="identity")+
  labs(y="Total reads",x="") +
  theme_pubr()+
  theme(legend.title=element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))
  

plot_phyla


 ggsave("Stacked_bar_total_abundance_site.pdf", width = 20, height = 15, units = "cm") 
#############
#STACKED BAR PLOTS 

#Merge by Sample type (column heading of your metadata). Don't worry about NAs.
merged_group_bh <-merge_samples(group_RA, "site")


otu_table(merged_group_bh)<-otu_table((otu_table(merged_group_bh))/(data.frame(sample_data(group_RA)) %>% group_by(site) %>% tally())$n,taxa_are_rows=FALSE)

#In this one you have the choice of leaving all taxa that are above x% in at least one sample in, or just merging everything low abundance on a sample by sample basis, keep_all=TRUE

#merged_phyla<-merge_low_abun_taxa(merged_phyla,1,keep_all=TRUE)

#Melt for use with ggplot
melted_group_bh<-psmelt(merged_group_bh)


plot_group_bh<-ggplot(data=melted_group_bh, 
                   aes(fill=Species,x=Sample, y=Abundance))+
  geom_bar(stat="identity")+
  labs(y="Number of reads",x="") +
  theme_pubr()+
  theme(legend.title=element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) 


plot_group_bh


ggsave(plot =plot_group_bh,"Stacked_bar_RA_depth.pdf", width = 20, height = 15, units = "cm") 


#FACATED bar plots from same melted data as above

facated_bar_tax <- ggplot(melted_phyla, aes(fill=Species, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity")+
  labs(y="Percent of Tetracladium reads",x="") +
  theme_pubr()+
  theme(legend.title=element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

facated_bar_tax
ggsave(plot =facated_bar_tax,"facated_bar_tax.pdf", width = 20, height = 15, units = "cm") 

############
# #SEM
# #Piecewise SEM ROOT
# min_max_norm <- function(x) {
#   (x - min(x)) / (max(x) - min(x))
# }
# #######data:
# #richness
# tryaa
# #reorder variables so cat variables are at the back
# mod_data <- tryaa[, c(1,2,6,7,8,9,10,11,12,13,14,15,16,3,4,5)]
# mod_data
# mod_data2<- na.omit(mod_data)
# 
# #apply Min-Max normalization
# mod_data_norm <- as.data.frame(lapply(mod_data2[1:13], min_max_norm))
# 
# #add cat variables back
# mod_data_norm$horizon <- mod_data2$horizon
# mod_data_norm$site <- mod_data2$site
# mod_data_norm$myco <- mod_data2$myco
# 
# mod_data_norm
# 
# 
# #RA
# Y
# Ysum <- rowSums(Y,na.rm = FALSE)
# 
# data_metaSEM <- data.frame(sample_data(phylo))
# mod_data_RA <- merge(Ysum, data_metaSEM, by ="row.names", all = TRUE )
# rownames(mod_data_RA) <- mod_data_RA[,1]
# mod_data_RA[,1] <- NULL
# 
# #reorder variables so cat variables are at the back
# mod_data_RA <- mod_data_RA[, c(1,2,6,7,8,9,10,11,12,13,14,15,16,3,4,5)]
# mod_data_RA
# mod_data_RA2<- na.omit(mod_data_RA)
# 
# #apply Min-Max normalization
# mod_data_RA_norm <- as.data.frame(lapply(mod_data_RA2[1:10], min_max_norm))
# 
# #add cat variables back
# mod_data_RA_norm$horizon <- mod_data_RA2$horizon
# mod_data_RA_norm$site <- mod_data_RA2$site
# mod_data_RA_norm$myco <- mod_data_RA2$myco
# 
# mod_data_RA_norm

# 
# Model_R <- psem(
#   lmer(Observed ~ ecm_tree.percent + ph +(1|site),mod_data_norm),
#   lmer(ecm_tree.percent~ ph + (1|site),mod_data_norm),
#   lm(c.percent ~ ecm_tree.percent, mod_data_norm))
# 
# SEM_R <- plot(
#   Model_R,
#   return = F,
#   node_attrs = data.frame(shape = "rectangle", style = "rounded", fixedsize = F, color = "#2F4F4F", fillcolor = "white", fontsize = 8, penwidth = 2),
#   edge_attrs = data.frame(style = "solid", color = "black", arrowhead = "none"),
#   ns_dashed = T,
#   alpha = 0.05,
#   show = "std",
#   digits = 3,
#   add_edge_label_spaces = F,
#   title = "Richness PiecewiseSEM final model")
###########
######Heatmap composition
#REMOVE taxa and samples with 0 no reads  
hm <- prune_taxa(taxa_sums(group_RA) > 0, group_RA)
hm <- prune_samples(sample_sums(hm)>0,hm)

p.famrel.heatmap <- plot_composition(hm,
                                     sample.sort = "neatmap",
                                     otu.sort = "neatmap",
                                     plot.type = "heatmap") +
  scale_fill_gradient(low = "white", high = "#E53F38", breaks = seq(min(hm@otu_table),max(hm@otu_table),0.1)) +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, 
                                   vjust = 0.85, 
                                   hjust = 1))

p.famrel.heatmap





ggsave(plot = p.famrel.heatmap, file = "heatmap.pdf", dpi = 300)
########## 
#network analyses
root_ps <- merge_phyloseq(AM, ECM, Tetracladium_ps)
sample_variables(root_ps)
#REMOVE taxa and samples with 0 no reads  
#root_ps <- prune_taxa(taxa_sums(root_ps) > 0, root_ps)
#root_ps <- prune_samples(sample_sums(root_ps)>0,root_ps)

filtered_obj <- conglomerate_taxa(root_ps, "Genus")
genus_co<- co_occurrence_network(filtered_obj, treatment = NULL,
                      classification = 'Genus')


genus_co
ggsave(plot = genus_co, file = "genus_co.pdf", scale = 2)

co_occurrence(filtered_obj, treatment = NULL, rho = 0.8, p = 0.05, cores = 1)


filtered_obj2 <- conglomerate_taxa(Tetracladium_ps, "Species")
genus_co2<- co_occurrence_network(filtered_obj2, treatment = NULL,
                                 classification = 'Species'  )


genus_co2
ggsave(plot = genus_co2, file = "genus_co.pdf", scale = 2)


AM_top = prune_taxa(names(sort(taxa_sums(AM_ra), TRUE))[1:20], AM_ra)
ECM_top = prune_taxa(names(sort(taxa_sums(ECM_ra), TRUE))[1:20], ECM_ra)
top_ps <- merge_phyloseq(AM_top, ECM_top, group_RA)

top_ps20 <- conglomerate_taxa(top_ps, "Genus")


top_co20<- co_occurrence_network(top_ps20, treatment = NULL,
                                  classification = 'Genus'  )

top_co20

co_net <- co_occurrence(top_ps20, treatment = NULL, rho = 0, p = 0.05, cores = 1)
co_net <- as.data.frame(co_net)
write_csv(co_net, "co_net.csv")

fdr1 <- p.fdr(
  pvalues = co_net$p,
  zvalues = "two.sided",
  threshold = 0.05,
  adjust.method = "BH",
  BY.corr = "positive",
  just.fdr = FALSE,
  default.odds = 1,
  estim.method = "set.pi0",
  set.pi0 = 1,
  hist.breaks = "scott",
  ties.method = "random",
  sort.results = FALSE,
  na.rm = TRUE)


plot(fdr1, ylim = c(0, 0.1))

ggsave(plot = top_co, file = "myco_top50_genus_co.pdf", scale = 1)

variable_correlation_network(group_RA, variables = c("ecm_tree.percent", "c.percent", "LitLCI", "ph", "Sand"),
                     treatment = NULL, subset = NULL,
                     classification = "Species", method = "pearson")

variable_correlation_heatmap(group_RA, variables= c("c.percent", "LitLCI", "ph", "Sand"),
                             treatment = "myco",
                             classification = 'Species', method = 'spearman', cores = 1,
                             colors = c("#2C7BB6", "white", "#D7191C"),
                             significance_color = 'black')


filtered_obj_all <- conglomerate_taxa(ps, "Genus")

common_taxa(filtered_obj_all, treatment = "site", subset = NULL, n = "all")

taxa_core_graph(group_RA, treatment = NULL, subset = NULL,
                frequencies = seq(0.1, 1, 0.1), abundance_thresholds = seq(0.01, 1, 0.01),
                colors = 'default',
                treatment_labels = NULL, sample_labels = "site", classification_labels= "Genus")

top = prune_taxa(names(sort(taxa_sums(ps_RA), TRUE))[1:50], ps_RA)
top_all <- merge_phyloseq(top, group_RA)

top_all2 <- conglomerate_taxa(top_all, "Genus")


top_co<- co_occurrence_network(top_all2, treatment = NULL,
                               classification = 'Genus')

top_co

co_net2 <- co_occurrence(top_all2, treatment = NULL, rho = 0, p = 0.05, cores = 1)
co_net2 <- as.data.frame(co_net2)
write_csv(co_net2, "co_net2.csv")


fdr2 <- p.fdr(
  pvalues = co_net2$p,
  zvalues = "two.sided",
  threshold = 0.05,
  adjust.method = "BH",
  BY.corr = "positive",
  just.fdr = FALSE,
  default.odds = 1,
  estim.method = "set.pi0",
  set.pi0 = 1,
  hist.breaks = "scott",
  ties.method = "random",
  sort.results = FALSE,
  na.rm = TRUE)


plot(fdr2, ylim = c(0, 0.1))

ggsave(plot = top_co, file = "all_top50_genus_co.pdf", scale = 1)



#########
#ENZIME
Enzime

Y_enz_ra <- veganotu(Enzime)
Ysum_enz_ra <- rowSums(Y_enz_ra,na.rm = FALSE)
dat_enzime <- data.frame(sample_data(Enzime))
mod_data_enz <- merge(Ysum_enz_ra, dat_enzime, by ="row.names", all = TRUE )
rownames(mod_data_enz) <- mod_data_enz[,1]
mod_data_enz[,1] <- NULL
mod_data_enz
mod_data_enz2 = subset(mod_data_enz, select = c(x,ecm_tree.percent, n.percent, c.percent,
                                                s.percent,ph,
                                                AP_Spring, AP_Summer, AP_Fall,
                                                NAG_Spring, NAG_Summer, NAG_Fall,
                                                BG_Spring,  BG_Summer,  BG_Fall,
                                                Phenox_Spring,   Phenox_Summer,   Phenox_Fall,
                                                Pheno_Spring, Pheno_Summer, Pheno_Fall))

mod_data_enz2<- na.omit(mod_data_enz2)
mod_data_enz2 <- mutate_all(mod_data_enz2, function(x) as.numeric(as.character(x)))

mod_data_enz_norm <- as.data.frame(lapply(mod_data_enz2, min_max_norm))

mod_data_enz_norm$CN <- mod_data_enz_norm$c.percent/mod_data_enz_norm$n.percent

mod_data_enz_norm <- mod_data_enz_norm %>% 
  filter_all(all_vars(!is.infinite(.)))


#stepwise variable selection
#  prevent fitting sub-models to different datasets
mod_enz<- lm(x~ecm_tree.percent+ 
             s.percent+ph+
             AP_Spring+ AP_Summer+ AP_Fall+
             NAG_Spring+ NAG_Summer+ NAG_Fall+
             BG_Spring+  BG_Summer+  BG_Fall+
             Phenox_Spring+   Phenox_Summer+   Phenox_Fall+
             Pheno_Spring+ Pheno_Summer+ Pheno_Fall +CN, mod_data_enz_norm)


options(na.action = "na.fail")
dd_enz<-dredge(mod_enz, trace = 2, extra = c("R^2", F = function(x)
  summary(x)$fstatistic[[1]]))

subset(dd_enz, delta < 4)

# Visualize the model selection table:
par(mar = c(3,5,6,4))
plot(dd_enz, labAsExpr = TRUE)
# Model average models with delta AICc < 4
model.avg(dd_enz, subset = delta < 4)

#or as a 95% confidence set:
model.avg(dd_enz, subset = cumsum(weight) <= .95) # get averaged coefficients

#'Best' model
summary(get.models(dd_enz, 1)[[1]])


best_mod_enz_ra <- lm(x ~ CN + ph  + 1, data = mod_data_enz_norm)

summary(best_mod_enz_ra)

CMplotR <- check_model(
  best_mod_enz_ra,
  dot_size = 2,
  line_size = 0.8,
  panel = TRUE,
  check = "all",
  alpha = 0.2,
  dot_alpha = 0.8,
  colors = c("#3aaf85", "#1b6ca8", "#cd201f"),
  theme = "see::theme_lucid",
  detrend = FALSE,
  verbose = TRUE)


#ENZIME RICHNESS
observed_D_enzime <- estimate_richness(Enzime_D, split = TRUE, measures = "Observed")
#observed_D_samples <- filter(observed_D_samples, Observed > 0)
tryaa_enzime<-cbind(observed_D_enzime,dat_enzime)

#select variables to keep (no NAs)
tryaa_enzime2 = subset(tryaa_enzime, select = c(Observed, ecm_tree.percent, n.percent, c.percent,
                                         s.percent,ph,
                                         AP_Spring, AP_Summer, AP_Fall,
                                         NAG_Spring, NAG_Summer, NAG_Fall,
                                         BG_Spring,  BG_Summer,  BG_Fall,
                                         Phenox_Spring,   Phenox_Summer,   Phenox_Fall,
                                         Pheno_Spring, Pheno_Summer, Pheno_Fall))


tryaa_norm<- na.omit(tryaa_enzime2)
tryaa_norm2 <- mutate_all(tryaa_norm, function(x) as.numeric(as.character(x)))

tryaa_norm2 <- as.data.frame(lapply(tryaa_norm2, min_max_norm))

# #add back site column
# tryaa_norm2<-tryaa_norm
# tryaa_norm2$site <- mod_data_cor$site

tryaa_norm2$CN <- tryaa_norm2$c.percent/tryaa_norm2$n.percent

tryaa_norm2 <- tryaa_norm2 %>% 
  filter_all(all_vars(!is.infinite(.)))

mod_enz_d<- lm(Observed~ecm_tree.percent+
               s.percent+ph+
               AP_Spring+ AP_Summer+ AP_Fall+
               NAG_Spring+ NAG_Summer+ NAG_Fall+
               BG_Spring+  BG_Summer+  BG_Fall+
               Phenox_Spring+   Phenox_Summer+   Phenox_Fall+
               Pheno_Spring+ Pheno_Summer+ Pheno_Fall + CN, tryaa_norm2)

dd_enz_d<-dredge(mod_enz_d, trace = 2, extra = c("R^2", F = function(x)
  summary(x)$fstatistic[[1]]))

summary(get.models(dd_enz_d, 1)[[1]])

best_mod_enz_d <- lm(Observed ~ CN + ph + 1, data = tryaa_norm2)
summary(best_mod_enz_d)

all_enz_plot <- plot_summs(best_mod_enz_ra, best_mod_enz_d,
                      scale = TRUE, robust = TRUE,
                      model.names = c("Observed","RA"))

ggsave("both_enz_sum.pdf", all_enz_plot, dpi = 300)

CMplotR <- check_model(
  best_mod_enz_d,
  dot_size = 2,
  line_size = 0.8,
  panel = TRUE,
  check = "all",
  alpha = 0.2,
  dot_alpha = 0.8,
  colors = c("#3aaf85", "#1b6ca8", "#cd201f"),
  theme = "see::theme_lucid",
  detrend = FALSE,
  verbose = TRUE)

####map
sites <- read_excel("~/Documents/FOREST_Ryan/MS1/sites.xlsx")

state <- map_data("state")

map1<- map_data("state") %>% 
  ggplot() +
  geom_polygon(aes(x=long, y=lat, group=group),colour = "gray85", fill = "gray80") +
  geom_point(data = sites, 
             aes(x = Long, y = Lat), 
             colour = "purple", alpha = 0.5) + 
  theme_minimal() +
  coord_fixed(1.3)

require(ggspatial)

map2<-map1 +  
  annotation_north_arrow(location = "br", 
                         which_north = "true", #north arrow
                         style = north_arrow_fancy_orienteering) 

ggsave("map1.pdf", map2, dpi = 300)

#number of reads

median_all = median(sample_sums(ps))
mean_all = mean(sample_sums(ps))
min_all =min(sample_sums(ps))
max_all = max(sample_sums(ps))
sd_all = sd(sample_sums(ps))

tetra_stat <- subset_taxa(ps, Genus=="g__Tetracladium")
min_tetra =min(sample_sums(tetra_stat))
max_tetra = max(sample_sums(tetra_stat))
sd_tetra = sd(sample_sums(tetra_stat))
