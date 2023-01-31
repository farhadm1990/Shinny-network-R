# Shinny-network
A function to create association network for microbiome data: bacteria-bacteria and feature-metabolite association.
This function is dedicated to make graph/network based on the spearman (also pearson) correlation and the significant level of this correlation corrected for false dicorevy rate (FDR) by Benjamini-Hochberg (by default, other methods are also accepted. See the help sheet for p.ajust() function). This is a costume function and as it doens't count for partial effects of taxa, you must only use it for visualization and not for validation of associations. The function is also able to perfomr these analysis with and without Centered-Log ratio (CLR) transformation to account for difference in read depth. For the input matrix, you can simply use the phyloseq object and the function will do the rest. By default, the graph will be made from a dataframe, based on the most significantly correlated ASVs. 

## Function parameters:
`clr`: if you want to centered-log ratio transform data. `rel_abund`: if you want to do cumulative normalisation on the data into their relative abundance. `Treatment_prune`: if you want to break your dataset into the sublevels of the treatment, default is "TRUE". `treatment`: the column in your dataset with which you want to breack your dataset into sub-groups of treatments. `treat_level`: the level of your interest in the treatment column, which you want to split the dataset for. `top_n_taxa`: specifies the number of most diverse taxa based on their standard deviation filtering_threshold: number of the samples each ASV should appear in order to be passed through the filter.
`cor_threshold`: is for filtering the taxa with the absolute value of correlation below the threshold value (0.55). `sig_threshold`: is the significant alpha for the adjusted p.value. `directed`: is used for making the edge atributes. The Default is FALSE.
`taxa_level`: is for making the network for different taxa, default is Genus. `color_pallet`: is the pallet from which the color of the vertice will be chosen. The dfault is `qual`, but it can also be 'div', 'qual', 'seq'. `edge_pos_color`: is the color given to the edges between ASVs with positive correlation. Default is cyan. `edge_neg_color`: is the color given to the edges between ASVs with negative correlation. Default is red. treatment_prune: if you want to prune your taxa according to a particular treatment. The defaule it FALSE. Therefore, the treatment argument will be ignored. treatment: if the treatment_prune argument is set to TRUE, treatment will be considered in the downstream analysis. the dfault is ct (control). But you can always cahnge it to one of your treatments. `y_cor`: by default is NULL, which means the association will made between taxa. If you, for instance, have a matrix of your metabolite data, e.g. SCFA, you can pass it onto it and the association will be build between taxa-taxa and taxa-metabolite.

```R

network_forger = function(data, rel_abund = TRUE, clr = TRUE, treatment_prune = FALSE, 
                          treatment = "treatment column name", treat_level = "you.name.it",
                          filtering_threshold = 30,  FDR_method = "BH", y_cor= NULL, 
                          cor_method = c("spearman", "pearson", "kendall"),
                          cor_threshold = 0.55, sig_threshold = 0.01,
                          directed = "FALSE", 
                          taxa_level = "Genus", top_taxa = 500, color_pallet = c("qual", "div", "seq"),
                          edge_pos_color = "cyan", edge_neg_color = "red"){

#We need to include the required packages into the function

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
for (pkg in c("devtools", "dada2", "phyloseq", "ALDEx2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}


if(!requireNamespace("pacman", quietly = TRUE)){
  install.packages("pacman")
    }

    library("pacman")
pacman::p_load(devtools, RColorBrewer,  phyloseq, tidyverse, igraph, visNetwork, microbiome, glue)

library(phyloseq)
library(igraph)
library(glue)
library(microbiome)


#Taking care of the data. Is it normalized, if not we must do it?
#Taking care of the data. Is it normalized, if not we must do it?
if(rel_abund == TRUE){

	if (sum(colSums(otu_table(data)))/ncol(otu_table(data)) == 100 ){
    #making a rel.abund df to be passed on vertices latter on
	data = data
	
     } else if (sum(colSums(otu_table(data)))/ncol(otu_table(data)) == 1) {

    data = data


    } else {

    data = transform_sample_counts(data, function(x) {x/sum(x)})

    }
	
} else if(rel_abund == FALSE) {

data = data

} else {

stop("Check your rel_abund entry")

}
    
#Now we must either do CLRT or not?

if(clr == TRUE){ #if we want to make the correlation matrix out of cetered-log ratio matrix
##Network based on cetnered log-ratio transformation
#There are two problems with building networkd with simple correlation approach:
#first, correlations can be distorted when applied to normalized or rarefied sequencing data and second, it is not clear how to select a
#meaningful threshold on their strength. In addition, the plotting code did not distinguish between positive and negative edges. As you know, normalization
#or rarefaction constrains the total sample sum and results in compositional data. The constraint can introduce correlations that are absent in the real data.
#However, we have to normalize or rarefy data, since otherwise correlations will be driven by sequencing depth differences.
#There are two ways to deal with this compositionality problem in microbial network construction.
#The first is to transform the data and the second is to use association measures that are not affected by compositionality. Here, we will look at the first option.
#In the CLR transform, each abundance value is divided by the geometric mean of its sample and then a logarithm is taken. The geometric mean is the Sth-root of the product
#of all values in a sample, where S is the number of the species. When you look at the CLR-transformed matrix, you will notice that it now contains negative values.
#These negative values rule out a number of association measures that assume data to be positive, such as the Bray-Curtis dissimilarity and the Kullback-Leibler dissimilarity.


geom.mean = function (x) {
    return (exp(mean(log(x))))
}
asv.tb = otu_table(data)
asv.clr = matrix(0, nrow = nrow(asv.tb), ncol = ncol(asv.tb))
for(i in 1:ncol(asv.clr)) {
    samps = asv.tb[,i] #first obtain the samples
    non.zeros = which(samps>0)
    g = geom.mean(samps[non.zeros])
    asv.clr [non.zeros, i] = log(samps[non.zeros]/g)
}
rownames(asv.clr) = rownames(asv.tb)
colnames(asv.clr ) = colnames(asv.tb)

otu_table(data ) <- otu_table(asv.clr, taxa_are_rows = TRUE)

} else if(clr == FALSE) {
data = data
} else {

stop("Check your clr-transformation entry")
}

#=============================================================#
#in higher taxonomic levels (e.g. species level) we have alot of duplicated names, e.g. "uncultured_bacterium" belonging to different genera. so we must make a unique name based on their gennus
#A function to create unique names for each ASV. It removes any NA in Order level then attempts to use the name of one level higher taxa for those
#who have similar names, e.g. uncultured_bacterium


# A custumized function for tax_glom()
#A function to create unique names for each ASV. It removes any NA in Order level then attempts to use the name of one level higher taxa for those
#who have similar names, e.g. uncultured_bacterium

gloomer = function(ps = data, taxa_level = taxa_level, NArm = "TRUE"){
    rank.names = c('Kingdom','Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')


#====================Sometimes in genus level, we might have multiple uncultured organisms, which if we want to make unique out of them for the species level it won't work====
    #since adding uncultured to uncultered is sill duplication. therefore if the taxa_level is set to species we first make a unique genus and then we go further to the speices===#

#Removing unculured Family
ps = subset_taxa(ps, !Family %in% c("uncultured", "NA"))

if(taxa_level == "Species") {
ps = subset_taxa(ps, !Genus %in% NA)#we remove genus tagged NA
phyloseq::tax_table(ps)[, taxa_level] <- ifelse(phyloseq::tax_table(ps)[, taxa_level] %in% NA, paste0("unknown"), paste(phyloseq::tax_table(ps)[, taxa_level]))#convert NA in species into unknown

   physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = NArm)
   taxdat = phyloseq::tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]

   taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
   otudat = otu_table(physeq)

#first take care of the uncultured genus
taxdat[,6] = ifelse(taxdat[,6] == "uncultured",
       paste0(taxdat[ , length(rank.names[1:which(rank.names=="Genus")])-1], "_", taxdat[,6]), paste(taxdat[,6]))

spec1 = taxdat[, taxa_level] %>% as.vector
spec2  = taxdat[, taxa_level] %>% as.vector

    uni  = matrix(NA, ncol = length(spec2), nrow = length(spec1))
    for(i in seq_along(spec1)){
    for(j in seq_along(spec2)){
    uni[i, j] = ifelse(spec1[i] == spec2[j] , "TRUE", "FALSE")
    }
        }

rownames(uni) <-spec1
colnames(uni) <- spec2
uni[upper.tri(uni, diag = TRUE)] = 0 #get rid of diagonals and upper triangle

duplis = uni %>% melt() %>% filter(value == "TRUE")

if(dim(duplis)[[1]] > 0) {
duplis = uni %>% melt() %>% filter(value == "TRUE") %>% dplyr::select(1) %>% unique() %>% unlist() %>% as.vector()
taxdat = taxdat %>% mutate( uni= ifelse(taxdat[, taxa_level] %in% duplis,
                    paste0(taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "_", taxdat[,taxa_level]), taxdat[,taxa_level]))

taxdat[, taxa_level] = taxdat[, "uni"]
taxdat[, "uni"] <- NULL
taxdat <- as.matrix(taxdat)
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- phyloseq::tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
phyloseq::tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat

} else {

taxdat <- as.matrix(taxdat)
taxdat <- phyloseq::tax_table(taxdat)
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- phyloseq::tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
phyloseq::tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat

}
#ps = phyloseq(otu_table(otudat, taxa_are_rows = T), tax_table(as.matrix(taxdat)), sample_data(physeq))



#==========================================#
} else if (taxa_level == "Genus") {

    physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = NArm)
    taxdat = phyloseq::tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]

   taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
    otudat = otu_table(physeq)

# take care of the uncultured genus
taxdat[,6] = ifelse(taxdat[,6] == "uncultured",
       paste(taxdat[ , length(rank.names[1:which(rank.names=="Genus")])-1], "_", taxdat[,6]), paste(taxdat[,6]))


rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[taxdat[,taxa_level] %in% rownames(otudat), taxa_level]
taxdat <- as.matrix(taxdat)
taxdat <- phyloseq::tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
phyloseq::tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
#ps = phyloseq(otu_table(otudat, taxa_are_rows = T), phyloseq::tax_table(as.matrix(taxdat)), sample_data(physeq))



} else {


physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = TRUE)
    taxdat = phyloseq::tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]

taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
otudat = otu_table(physeq)

spec1 = taxdat[, taxa_level] %>% as.vector
spec2  = taxdat[, taxa_level] %>% as.vector

    uni  = matrix(NA, ncol = length(spec2), nrow = length(spec1))
    for(i in seq_along(spec1)){
        for(j in seq_along(spec2)){
    uni[i, j] = ifelse(spec1[i] == spec2[j] , "TRUE", "FALSE")
    }
        }

rownames(uni) <-spec1
colnames(uni) <- spec2
uni[upper.tri(uni, diag = TRUE)] = 0 #get rid of diagonals and upper triangle

duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE")

if(dim(duplis)[[1]] > 0){#if there is not duplications, we can simply use the taxa names as the row name

    duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE") %>% dplyr::select(1)%>% unique() %>% unlist %>% as.vector
taxdat = taxdat %>% mutate( uni= ifelse(taxdat[, taxa_level] %in% duplis,
                    paste(taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "_", taxdat[,taxa_level]), taxdat[,taxa_level]))

taxdat[, taxa_level] = taxdat[, "uni"]
taxdat[, "uni"] <- NULL
taxdat <- as.matrix(taxdat)
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- phyloseq::tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
phyloseq::tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
} else {

taxdat <- as.matrix(taxdat)
taxdat <- phyloseq::tax_table(taxdat)
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- phyloseq::tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
phyloseq::tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
}
#ps = phyloseq(otu_table(otudat, taxa_are_rows = T), tax_table(as.matrix(taxdat)), sample_data(physeq))


}
return(physeq)
    }
#=================================================================================================#
#=================================================================================================#

#Aglomerating the taxa to the specified taxonomic level
physeq = gloomer(data, taxa_level = taxa_level, NArm = TRUE)

#sometimes you want to make the graph for only one level of your treatments, e.g. only for healthy gourps or only for gorups treated with DSS
if(treatment_prune == FALSE){
   
	asv = otu_table(physeq)


} else if(treatment_prune == "TRUE"){


    vect.id = sample_data(physeq)[,treatment] %>% pull %>% unfactor == treat_level
    physeq = prune_samples(vect.id, physeq)
	
    asv = otu_table(physeq)
} else {
    
stop("Your treatment_prune entry must be TRUE/FALSE") 
  
}


#3. Building the correlation matrix


asv = asv[!rownames(asv) %in% c("NA", "Unknown","", "uncharacterized", "unassigned", "unknown", "uncultured"), ] #to remove Unknown bacteria
asv.relabund = asv %>% data.frame%>% mutate(rel.abund  = rowSums(asv)/sum(asv)*100) 

#Filtering the ASVs who appeared in less than 30% of the samples out of the dataset
threshold = filtering_threshold
index = asv
index[asv > 0] = 1

asv.filt = asv[rowSums(index)/ncol(index)*100 >= threshold, ]

greacer <- function(ps, comp, treat, method , top_n_taxa,  y, pseudo = -1000000, FDR){



    
treat.vect <- sample_data(ps)[, treat] %>% pull %>% unique %>% as.vector

    
#Making a correlation matrix
if(is.null(y)){
    
  dat.ls <- list()
  for(i in treat.vect){
      
       dat.ls[[i]] <- comp[, colnames(comp) %in% rownames(ps@sam_data[ps@sam_data[,treat]== i])] 
       
      
}
    cor_main <- list()
    
    for(i in treat.vect){
                
        cor_main[[i]] <- Hmisc::rcorr(x = t(dat.ls[[i]] ), y = t(dat.ls[[i]]), type = method)
    
}

cor.pval <- list()
Cor <- list()
 sds <- list()  
    ord <- list()
    cor.sds <- list()  
    ord.cor <- list()
for(i in treat.vect){
   cor.pval[[i]] = cor_main[[i]]$P
    cor.pval[[i]][upper.tri(cor.pval[[i]], diag = TRUE)] <- pseudo
    
    
    Cor[[i]] = cor_main[[i]]$r 
    cor.sds[[i]] <- rowSds(Cor[[i]], na.rm = TRUE)
    ord.cor[[i]] <- order(cor.sds[[i]], decreasing = TRUE)[1:top_n_taxa]
    Cor[[i]] <- Cor[[i]][ord.cor[[i]], ord.cor[[i]]]
    Cor[[i]][upper.tri(Cor[[i]], diag = TRUE)] <- pseudo
    
    cor.pval[[i]] = cor.pval[[i]][rownames(cor.pval[[i]]) %in% rownames(Cor[[i]]), colnames(cor.pval[[i]]) %in% colnames(Cor[[i]])]
    }

cor.pval = cor.pval %>%  
   reshape2::melt(value.name = "pval", 
   varnames = c("from", "to")) %>% filter(pval != pseudo)
cor.pval = cor.pval[complete.cases(cor.pval),]
    
Cors = Cor%>% 
     reshape2::melt(value.name = "cor", 
      varnames = c("from", "to")) %>% filter(cor != pseudo)
Cors = Cors[complete.cases(Cors),]
        
    } else {
    
 
    
  if (identical(sort(colnames(y)), sort(colnames(comp)))){
    y = as.matrix(y)
} else {
    y = t(as.matrix(y))
}
#check if they have the same dimensions
    if(dim(comp)[2] > dim(y)[2]){
    
    comp = comp[, colnames(comp)%in% colnames(y)]
    
} else if (dim(comp)[2] < dim(y)[2]){
    y = y[, colnames(y)%in% colnames(comp)]
} else {
    comp = comp
    y = y
}
    
   dat.ls <- list()
    dat.ls.y <- list()
    
  for(i in treat.vect){
      
       dat.ls[[i]] <- comp[, colnames(comp) %in% rownames(ps@sam_data[ps@sam_data[,treat]== i])]  
       dat.ls.y[[i]] <- y[, colnames(y) %in% rownames(ps@sam_data[ps@sam_data[,treat]== i])]
      
}
    cor_main <- list()
    
    for(i in treat.vect){
                
        cor_main[[i]] <- Hmisc::rcorr(x = t(dat.ls[[i]] ), y = t(dat.ls.y[[i]]), type = method)
        
    }
   
cor.pval <- list()
Cor <- list() 
  #sds <- list()
  #  ord <- list()
   # cor.sds <- list()
    #cor.ord <- list()
for(i in treat.vect){
    cor.pval[[i]] = cor_main[[i]]$P[!rownames(cor_main[[i]]$P) %in% rownames(dat.ls.y[[i]]),] 
    #sds[[i]] <- rowSds(cor.pval[[i]], na.rm = TRUE)
    #ord[[i]] <- order(sds[[i]], decreasing = TRUE)[1:top_n_taxa]
    #cor.pval[[i]] <- cor.pval[[i]][ord[[i]],]
    
    Cor[[i]] = cor_main[[i]]$r[!rownames(cor_main[[i]]$r) %in% rownames(dat.ls.y[[i]]),]
   # cor.sds[[i]] <- rowSds(Cor[[i]], na.rm = TRUE)
    #cor.ord[[i]] <- order(Cor[[i]], decreasing = TRUE)[1:top_n_taxa]
   # Cor[[i]] <- Cor[[i]][cor.ord[[i]],]
}
 cor.pval = cor.pval %>%  reshape2::melt(value.name = "pval", varnames = c("from", "to"))
    #cor.pval = cor.pval[complete.cases(cor.pval),]
 Cors = Cor %>%  reshape2::melt(value.name = "cor", varnames = c("from", "to"))
    #Cors = Cors[complete.cases(Cors),]
}
    
  
   cor.pval[, dim(cor.pval)[2]]<- factor(cor.pval[,dim(cor.pval)[2]])
   colnames(cor.pval)[dim(cor.pval)[2]] <- treat    
   cor.pval$q.vals <- p.adjust(cor.pval$pval, method = FDR)
    
   Cors[, dim(Cors)[2]]<- factor(Cors[,dim(Cors)[2]])
   colnames(Cors)[dim(Cors)[2]] <- treat

  long.cor = Cors
  long.pval = cor.pval

long.cor$index <- paste(long.cor$from, long.cor$to)
long.pval$index <- paste(long.pval$from, long.pval$to)

combined.df <- left_join(long.cor, long.pval %>% 
          dplyr::select(q.vals, index), by = "index")%>% group_by(index) %>% 
		  distinct(.[, treat] , .keep_all = T)  %>% arrange(desc(index))
combined.df = combined.df[complete.cases(combined.df),]
    
combined.df$index <- NULL 
    

    
return(combined.df)    
}



treat_level =treatment 
comp = asv.filt
FDR = FDR_method
ps = physeq
top_taxa =top_taxa
method = cor_method
y = y_cor
#4. making the graph from dataframe

cor.filt = greacer(y = y, top_n_taxa = top_taxa, ps = ps, method = method, comp = comp, treat = treat_level, FDR = FDR)

edge =cor.filt[abs(cor.filt$cor)>cor_threshold,]#we only choose the strongest correlations (> 0.55)
edge = cor.filt[cor.filt$q.vals <= sig_threshold, ]
edge$from <-edge$from %>% unfactor %>% factor #to relevel to the new levels
edge$to <- edge$to %>% unfactor %>% factor

#Sanity check to remove correlation of ASVs to themselves
edge = edge[!unfactor(edge$from) == unfactor(edge$to),]


#making the network 
graph.df = graph_from_data_frame(edge, directed = FALSE)

#removing the multiple edges from the graph, i.e. those who go twice
graph.df = delete_edges(graph = graph.df, edges = E(graph.df)[which_multiple(graph.df)])

#Vertice attributes
V(graph.df)$degree =igraph::degree(graph.df, mode = "all")
V(graph.df)$label = V(graph.df)$name #create label for the network
V(graph.df)$eigen <- igraph::evcent(graph.df)$vector
V(graph.df)$betweeness <- betweenness(graph.df, directed=FALSE)
V(graph.df)$rel.abund <- asv.relabund$rel.abund[rownames(asv.relabund) %in% names(V(graph.df))]
#V(graph.df)$treat <- levels(node_met$treatment)
#V(graph.df)$sampleID <- rownames(node_met)

#Edge atributes
E(graph.df, directed = directed)$direction <- ifelse(E(graph.df)$cor>0, "Pisitive", "Negative")
E(graph.df, directed = directed)$color <- ifelse(E(graph.df)$cor>0, edge_pos_color, edge_neg_color)
E(graph.df, directed = directed)$weight =round(abs(E(graph.df)$cor), 2)

edge$cor <- round(edge$cor, 2)

#making a name index for the phylums to be painted on the vectors
taxa.names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
color.index = c("King.col", "Phyl.col", "Class.col", "Ord.col", "Fam.col", "Gen.col", "Spec.col")
    
#making a color vector for the taxa to be painted on the vertices
library(RColorBrewer)

#first make a matrix with colors and taxa level
pals = brewer.pal.info[brewer.pal.info$category %in% color_pallet,]
col_vector = unlist(mapply(brewer.pal, pals$maxcolors, rownames(pals))) %>% unique

if(taxa_level == "Species"){
    
    taxa.table = taxdat[complete.cases(taxdat), seq(which(taxa.names %in% taxa_level))]
    

} else { 
    taxa.table = tax_table(physeq)[complete.cases(tax_table(physeq)),seq(which(taxa.names %in%taxa_level))]
    taxa.table = taxa.table[!rownames(taxa.table)%in%"Unknown",]
}
    
    
colindex=matrix(NA, nrow = nrow(taxa.table), ncol = ncol(taxa.table))#making different taxa dataframe correspondent to the name of our "From" vertex names.
 for(j in 1:ncol(colindex)){
    colindex[,j] = col_vector[as.numeric(as.factor(taxa.table[, j]))] 
      
 } 

colnames(colindex) = color.index[seq(ncol(colindex))]
rownames(colindex) = rownames(taxa.table)
colindex = cbind(taxa.table, colindex)

if(taxa_level == "Kingdom"){
    

for(n in taxa.names[1:ncol(taxa.table)]){
            for(m in color.index[1:ncol(taxa.table)]){
                
            vertex_attr(graph.df, n)= colindex[rownames(colindex) %in% names(V(graph.df)),n] %>% as.vector
            vertex_attr(graph.df, m)= colindex[rownames(colindex) %in% names(V(graph.df)),m] %>% as.vector
            
 
}
    } 
    
} else if (taxa_level == "Phylum") {

for(n in taxa.names[1:ncol(taxa.table)]){
            for(m in color.index[1:ncol(taxa.table)]){
                
            vertex_attr(graph.df, n)= colindex[rownames(colindex) %in% names(V(graph.df)),n] %>% as.vector
            vertex_attr(graph.df, m)= colindex[rownames(colindex) %in% names(V(graph.df)),m] %>% as.vector
            
 
}
    } 
  
} else if(taxa_level == "Order"){
    
 for(n in taxa.names[1:ncol(taxa.table)]){
            for(m in color.index[1:ncol(taxa.table)]){
                
            vertex_attr(graph.df, n)= colindex[rownames(colindex) %in% names(V(graph.df)),n] %>% as.vector
            vertex_attr(graph.df, m)= colindex[rownames(colindex) %in% names(V(graph.df)),m] %>% as.vector
            
 
}
    }  
    
} else if(taxa_level == "Class") {
  
    for(n in taxa.names[1:ncol(taxa.table)]){
            for(m in color.index[1:ncol(taxa.table)]){
                
            vertex_attr(graph.df, n)= colindex[rownames(colindex) %in% names(V(graph.df)),n] %>% as.vector
            vertex_attr(graph.df, m)= colindex[rownames(colindex) %in% names(V(graph.df)),m] %>% as.vector
            
 
}
    } 
    
} else if(taxa_level == "Family") {
 
    for(n in taxa.names[1:ncol(taxa.table)]){
            for(m in color.index[1:ncol(taxa.table)]){
                
            vertex_attr(graph.df, n)= colindex[rownames(colindex) %in% names(V(graph.df)),n] %>% as.vector
            vertex_attr(graph.df, m)= colindex[rownames(colindex) %in% names(V(graph.df)),m] %>% as.vector
            
 
}
    } 
    
} else if(taxa_level == "Genus"){

for(n in taxa.names[1:ncol(taxa.table)]){
            for(m in color.index[1:ncol(taxa.table)]){
                
            vertex_attr(graph.df, n)= colindex[rownames(colindex) %in% names(V(graph.df)),n] %>% as.vector
            vertex_attr(graph.df, m)= colindex[rownames(colindex) %in% names(V(graph.df)),m] %>% as.vector
            
 
}
    } 
     
} else if(taxa_level == "Species"){
    

    for(n in taxa.names[1:ncol(taxa.table)]){
            for(m in color.index[1:ncol(taxa.table)]){
                
            vertex_attr(graph.df, n)= colindex[rownames(colindex) %in% names(V(graph.df)),n] %>% as.vector
            vertex_attr(graph.df, m)= colindex[rownames(colindex) %in% names(V(graph.df)),m] %>% as.vector
            
 
}
    } 


 } else {
print('Warning! the taxa name has not been entered correctly :( 
\n it must be one of this list: "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"')
}

 #since the labels of the vertices are automatically oredered alphabetically and this would make problem with assigning right color to them,
    #we make the labels exactly similar to the taxa_level
v.labs = factor(vertex_attr(graph.df, taxa_level), levels =  vertex_attr(graph.df, taxa_level))

graph.df = igraph::set.vertex.attribute(graph.df, "label", value=paste(v.labs %>% levels))

print(glue::glue("Your graph is ready! Graph is based on {cor_method} correlation matrix\n for CLR-transformed data, chosen based on BH q.value <= {sig_threshold}")     )



structure(list(asv.filt = asv.filt, cor.pval = cor.pval, cor = cor, edge.df = edge,
                   graph = graph.df , cor.filt = cor.filt, colors = colindex), class = "GraphPack")





}




```

### Testing the function

```R
graph.dig.nodiar = network_forger(data = gen.dig, sig_threshold = 0.01, cor_threshold = 0.6, rel_abund = F, 
				  clr = T, taxa_level = "Genus", y_cor = NULL, cor_method = "spearman",
                                  top_taxa = 100, treatment_prune = TRUE, treat_level = "NoDiar", treatment = "status")

```
### Visulizing the network
```R
ggnetwork::ggnetwork(graph.dig.nodiar$graph, layout = layout.kamada.kawai(graph.dig.nodiar$graph)) %>% 
ggplot(aes(x = x, y = y, xend = xend, yend = yend)) + 
theme_blank()+ 
    geom_edges(aes(color = direction)) + 
geom_nodes(aes(fill = Phylum), shape = 21, size = 10) + 
geom_nodetext_repel(aes(label = Genus), nudge_x = 0, 
                    nudge_y = 0.05, size = 5,point.padding = unit(0.5, "lines"))+
#geom_edgetext(aes(label = round(cor, 2)), size = 1, color = "black", fill = "white", label.size = 0) +
labs(color = "Correlation") + scale_color_manual(values = c("red", "cyan")) + 
ggtitle("NoDiar, cor > 0.6, FDR <= 0.1")
```

![Association network of taxa](https://github.com/farhadm1990/Shinny-network/blob/main/NoDiar.jpeg)
> Figure 1. Significant (p.adjust < 0.01) association network between different Genera based on Spearman rank test. Color of each vertex is the associated phyla and the edge color is positive or negative correlation. 

### Network of taxa with chemical biomarkers
```R
comp =asv.filt
y  = cor.scfa   
 
    
if (identical(sort(colnames(y)), sort(colnames(comp)))){
    y = as.matrix(y)
} else {
    y = t(as.matrix(y))
}
#check if they have the same dimensions
    if(dim(comp)[2] > dim(y)[2]){
    
    comp = comp[, colnames(comp)%in% colnames(y)]
    
} else if (dim(comp)[2] < dim(y)[2]){
    y = y[, colnames(y)%in% colnames(comp)]
} else {
    comp = comp
    y = y
}
    
   dat.ls <- list()
    dat.ls.y <- list()
    
  for(i in treat.vect){
      
       dat.ls[[i]] <- comp[, colnames(comp) %in% rownames(ps@sam_data[ps@sam_data[,treat]== i])]  
       dat.ls.y[[i]] <- y[, colnames(y) %in% rownames(ps@sam_data[ps@sam_data[,treat]== i])]
      
}
    cor_main <- list()
    
    for(i in treat.vect){
                
        cor_main[[i]] <- Hmisc::rcorr(x = t(dat.ls[[i]] ), y = t(dat.ls.y[[i]]), type = method)
        
    }
   
cor.pval <- list()
Cor <- list() 
  #sds <- list()
  #  ord <- list()
   # cor.sds <- list()
    #cor.ord <- list()
for(i in treat.vect){
    cor.pval[[i]] = cor_main[[i]]$P[!rownames(cor_main[[i]]$P) %in% rownames(dat.ls.y[[i]]),] 
    #sds[[i]] <- rowSds(cor.pval[[i]], na.rm = TRUE)
    #ord[[i]] <- order(sds[[i]], decreasing = TRUE)[1:top_n_taxa]
    #cor.pval[[i]] <- cor.pval[[i]][ord[[i]],]
    
    Cor[[i]] = cor_main[[i]]$r[!rownames(cor_main[[i]]$r) %in% rownames(dat.ls.y[[i]]),]
   # cor.sds[[i]] <- rowSds(Cor[[i]], na.rm = TRUE)
    #cor.ord[[i]] <- order(Cor[[i]], decreasing = TRUE)[1:top_n_taxa]
   # Cor[[i]] <- Cor[[i]][cor.ord[[i]],]
}
 cor.pval = cor.pval %>%  reshape2::melt(value.name = "pval", varnames = c("from", "to"))
    #cor.pval = cor.pval[complete.cases(cor.pval),]
 Cors = Cor %>%  reshape2::melt(value.name = "cor", varnames = c("from", "to"))
    #Cors = Cors[complete.cases(Cors),]
}
```

![Association graph between taxa and chemicals](https://github.com/farhadm1990/Shinny-network-R/blob/main/Graph_taxa_biomarker.png)
> Figure 2. Significant (p.adjust < 0.01) association network between different Genera based on Spearman rank test and chemicals produced by bacteria. Color of each vertex is the associated treatment and the edge color is positive or negative correlation. 
