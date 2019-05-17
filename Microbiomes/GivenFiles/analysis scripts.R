#install Phyloseq

library("phyloseq")
getwd()
#Import data Fill=TRUE reads taxonomy with some missing values
OTU=read.table("~/Dropbox/Keck/cmaiki_plot1/Amend_run_Dec7/abundance_table_annotated_ID=100.tsv", row.names=1, header=TRUE, check.names=FALSE, fill=TRUE)
OTU=OTU[,1:(ncol(OTU)-1)]#delete taxonomy column
names(OTU) = sub("WMEA_","1",names(OTU)) #delete prefix add a 1
names(OTU) = sub("_.*","",names(OTU))  #delete postfix
phyloseq_OTU=otu_table(OTU, taxa_are_rows = TRUE)
#had to separate taxonomies manualy in excel
ID100_taxonomy <- as.matrix(read.delim("~/Dropbox/Keck/cmaiki_plot1/Amend_run_Dec7/constax_otu100/outputs/consensus_taxonomy.txt", row.names=1, stringsAsFactors=FALSE))
TAX = tax_table(ID100_taxonomy)
PSsample <- read.csv("~/Dropbox/Keck/cmaiki_plot1/CMAIKI_Miseq1MappingFile_20181105.csv", row.names=1, header=TRUE,stringsAsFactors=FALSE)
sampledata=sample_data(PSsample)
physeq1=merge_phyloseq(sampledata,TAX,phyloseq_OTU)
#######
#######
#######
#######
library("ggplot2")
library("plyr")
library("vegan")
theme_set(theme_bw())
rank_names(physeq1)
sample_sums(physeq1)
PS_1000=prune_samples(sample_sums(physeq1)>1000, physeq1)
PS_frac=transform_sample_counts(PS_1000, function(OTU) OTU/sum(OTU) )
plot1_ord=ordinate(PS_frac,"NMDS","bray")
p1=plot_ordination(PS_frac, plot1_ord,type="samples",color="taq", shape="trophic")
print(p1)
##########
##########
##########
DS=subset_samples(physeq1, collection_label=="Drosophila")
plot_bar(DS,fill="Order")
PS=filter_taxa(PS, function(x) sum(x > 3) > (0.1*length(x)), TRUE)#filter taxa seen fewer than 3 times in 1% of the samples

PS1000=transform_sample_counts(PS, function(OTU) OTU/sum(OTU) )
plot_heatmap(PS1000, "NMDS", "bray", sample.label = "sample_type", taxa.label="F")
heatmap(otu_table(PS1000))

##########
##########
##########

Merged_sample_type=merge_samples(physeq1, "sample_type") 
Merged_sample_type
Merged_troph=merge_samples(physeq1, "trophic") 
Merged_troph
troph_otu=as(otu_table(Merged_troph), "matrix")
temptroph=nestedtemp(troph_otu)
temptroph
plot(temptroph)  
sample_otu=as(otu_table(Merged_sample_type), "matrix")
sample_otu=nestedtemp(troph_otu)
plot(sample_otu)
plot_richness(Merged_troph)
rar=rarefy_even_depth(physeq1, 2000)
plot_richness(rar, x="sample_type", color="taq",measures=c("Observed", "Shannon"))
plot_richness(Merged_sample_type, sortby="Observed",measures=c("Observed", "Shannon"))
?plot_richness


############
############
#############
families=tax_glom(physeq1, taxrank=rank_names(physeq1)[5],NArm=TRUE)
families=prune_samples(sample_sums(families)>500, families)
famnorm=transform_sample_counts(families, function(OTU) OTU/sum(OTU) )
plot_heatmap(families, "NMDS", "bray", sample.label = "sample_type", taxa.label="Family", sample.order="sample_type")
