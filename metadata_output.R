---
  title: "OTU Workflow"
author: "Sarah Christofides"
date: "25/10/2021"
output: html_document
---
  
  For a general introduction to statistics and working in R, go to http://environmentalcomputing.net/ - **highly recommended**! A useful (although more technical) guide to multivariate statistics in particular is found at https://sites.google.com/site/mb3gustame/ (GUSTA ME). 

It is recommended that you save and edit a separate copy of this script for every dataset you analyse. Always keep a 'clean' master copy. Don't have just one 'working' copy which you alter for each subsequent dataset - create  a new working copy for each new set of data. This will ensure that your past analyses remain reproducible.

Some of the techniques are demonstrated for just one predictor variable. Obviously, you may wish to apply them to more than one predictor. In this case copy, paste and edit chunks of code as needed. As for the script as a whole, DON'T just overwrite bits of code/use the same object names for different things. (For example, if you wish to run a second PERMANOVA call it *pma2* or similar rather than calling both *pma*.) This is a recipe for disaster!!!
  
  ## Set up the session
  
  ### Set the working directory
  
  ```{r}
setwd("~/mydata/proj/MiniProject")
```
### Install packages

I have included this for the sake of completion, but in fact all the packages you need should be pre-installed on the server. For this reason I have commented out the commands below.

```{r}
devtools::install_github("jbisanz/qiime2R")
install.packages("metacoder")
install.packages("randomForest")
install.packages("vegan")

#Install metacoder development package

install.packages("devtools")



download.file("https://cran.r-project.org/src/contrib/Archive/bold/bold_1.3.0.tar.gz", destfile = "bold_1.3.0.tar.gz")

install.packages("bold_1.3.0.tar.gz", type="source", repos=NULL)

remotes::install_github("KlausVigo/phangorn")
remotes::install_github("ropensci/taxize")
remotes::install_github("crul/bold")
devtools::install_github("grunwaldlab/metacoder")



library(metacoder)
```

### Load packages (you need to do this every time you open R and want to use them)

```{r}
library(qiime2R)
library(metacoder)
library(randomForest)
library(vegan)
```

### Set up the colours you're going to use - run this once without changing anything

This is a palette of ten colours selected for being colour-blind friendly.

```{r}
cbfcol<-c(
  light.orange=rgb(red=255,green=191,blue=127,maxColorValue=255),
  bright.orange=rgb(red=255,green=127,blue=0,maxColorValue=255),
  light.green=rgb(red=178,green=255,blue=140,maxColorValue=255),
  bright.green=rgb(red=50,green=255,blue=0,maxColorValue=255),
  light.blue=rgb(red=165,green=237,blue=255,maxColorValue=255),
  bright.blue=rgb(red=25,green=178,blue=255,maxColorValue=255),
  lilac=rgb(red=204,green=191,blue=255,maxColorValue=255),
  purple=rgb(red=101,green=76,blue=255,maxColorValue=255),
  pink=rgb(red=255,green=153,blue=191,maxColorValue=255),
  red=rgb(red=229,green=25,blue=50,maxColorValue=255),
  black=rgb(red=0,green=0,blue=0,maxColorValue=255)
)
palette(cbfcol)

```

# Prepare the data for analysis

## Bacterial data preparation

We will use a very useful package called qiime2R to read the ASVs in directly from the qza. We will also need a table of metadata.

```{r}
wholedata<-read_qza("table.qza")

#Separate just the bacteria (ASV) columns
bact<-wholedata$data
#Make the ASVs rows and the samples columns 
bact<-as.data.frame(t(bact))
#Remove the wholedata object to free up memory space
rm(wholedata)

#Read in the metadata columns
metadata<-read.delim("sample-metadata2.tsv", stringsAsFactors = T)

#VERY IMPORTANT!! Check the the ASV table and metadata have samples in the same order.
all.equal(rownames(bact), as.character(metadata$sample.id))
#Put the data in the same order as the metadata
bact<-bact[match(metadata$sample.id, rownames(bact)),]
#Check again
all.equal(rownames(bact), as.character(metadata$sample.id))

```

### Check and clean your metadata

It is important to check that your metadata file is in good shape

```{r}
summary(metadata)

#CHECK: are there any factors which have been read in as numeric? (This will show up in the summary by having min/mean/max rather than the number of observations for each factor level). Example:
metadata$year<-as.factor(metadata$year)

#CHECK: are there any misspellings/ inconsistent capital letters in the factors? Example:
metadata$body.site[metadata$body.site=="Gut"]<-"gut"
metadata$body.site<-droplevels(as.factor(metadata$body.site))

#CHECK: are the factors in the order you would like them to be? Example: in a yes/no factor 'no' would be default be the reference level, as it comes first alphabetically. To make 'yes' the reference level:
metadata$reported.antibiotic.usage<-relevel(metadata$reported.antibiotic.usage, ref = "Yes")
#Change the order of the levels in body site
metadata$body.site<-factor(metadata$body.site, levels = c("left palm", "right palm", "gut", "tongue"))
```

### Calculate the sequencing depth for bacterial samples

NGS sequencing produces uneven numbers of reads for each sample. We are going to explore how many reads each sample has.

```{r}
metadata$depth<-rowSums(bact)

#Look at the range of sequence depth
metadata$depth[order(metadata$depth)]
#See which (if any) samples have less than 2000 reads
metadata[metadata$depth < 2000,]

#Remove samples below threshold of 2000 reads
#bact2<-bact[metadata$depth > 2000,]
#metadata2<-metadata[metadata$depth > 2000,]

#Convert depth to a proportion of the maximum for later use in plotting
metadata$grad<-metadata$depth/max(metadata$depth)
```

## Explore the bacterial data

### Plot a rarefaction curve

Low step number = detailed but slow, high step number = coarse but fast

```{r}
rarecurve(bact, step=250, col=cbfcol[as.numeric(metadata$body.site)], lwd=2, ylab="OTUs", label=FALSE)
legend("topright", inset=c(0,0.05),  levels(metadata$body.site), lwd=3, bty="n", cex=0.8, col=cbfcol)
#dev.print(pdf, "rarefaction_curve.pdf", height=6,width=9)
```

# Bray-Curtis NMDS

NMDS is an ordination technique - a way to visualise highly-dimensional data in a lower-dimensional way. NMDS requires use of a distance metric, and Bray-Curtis is commonly used for community data.

The quality of the fit is assessed by the stress: <0.05 is excellent, 0.05-0.1 is very good, 0.1-0.2 is adequate, 0.2 - 0.3 can cause problems, > 0.3 you need to increase the number of dimensions (via 'k=') that you allow the ordination to be carried out in.

```{r}
BCnmds<- metaMDS(bact, trymax = 500, k=2, pc="TRUE", 
                 distance = "bray", autotransform = FALSE)
BCnmds


#Plot coloured by sequencing depth
plot(BCnmds$points[,2]~BCnmds$points[,1], pch=20, col=rgb(red=1,blue=0,green=0,alpha=metadata$grad), xlab="NMDS1", ylab="NMDS2")
#dev.print(pdf, "NMDS_depth.pdf")
```

## Plot by each of the predictors

```{r}
#Colour by body site
plot(BCnmds$points[,2]~BCnmds$points[,1], pch=20, col=cbfcol[as.numeric(metadata$body.site)], xlab="NMDS1", ylab="NMDS3")
legend("bottomleft",  levels(metadata$body.site), pch=c(20), bty="n", cex=0.8, col=cbfcol)
ordiellipse(BCnmds$points, metadata$body.site, choices=c(1,2), conf=0.95, col=cbfcol)

#Try expanding the axes
plot(BCnmds$points[,2]~BCnmds$points[,1], pch=20, col=cbfcol[as.numeric(metadata$body.site)], xlab="NMDS1", ylab="NMDS3", xlim=c(-6, 4), ylim=c(-2,3))
legend("bottomleft",  levels(metadata$body.site), pch=c(20), bty="n", cex=0.8,
       col=cbfcol)
ordiellipse(BCnmds$points, metadata$body.site, choices=c(1,2), conf=0.95, col=cbfcol)

#Add in other predictors you wish to plot!
```

# PERMANOVA

PERMANOVA is the technique you'll use to do formal hypothesis testing on your multivariate dataset, as you would use ANOVA or linear regression on univariate data. The inner workings are different, but the way you specify the input and interpret the output is essentially the same. Like a GLM, PERMANOVA can handle a mix of categorical and continuous predictor (independent/ x) variables. If you aren't sure what those terms mean, now is the time to look them up!
  When doing ANOVA *etc.*, you would check for normality and equal variances. In contrast, PERMANOVA has only one formal assumption: that the disperson (spread of points) is not signifiantly different between groups.
In this example code, I assume you have two predictors: stream and surface. I'll assume you want to test for an effect of each one and their interaction. Obviously for your own data you will need to decide what hypotheses *you* want to test.

This example uses *method="bray"* which is suitable for community data.

### Check dispersion is equal between groups

Whereas for an ordinary ANOVA you would check for normality and homoscedasticity, PERMANOVA has only one formal assumption: that the dispersions within treatment groups are not significantly different.

```{r}
datdist<-vegdist((bact), method="bray")
#Check that the site groups AREN'T significantly different to each other
TukeyHSD(betadisper (datdist,metadata$body.site))
#Check that the time groups AREN'T significantly different to each other
TukeyHSD(betadisper (datdist,metadata$year))
#Check that the subject groups AREN'T significantly different to each other
TukeyHSD(betadisper (datdist,metadata$subject))
```

If the dispersions are not equal between groups, try using a fourth root transformation to squash abundance differences.

```{r}
froot<-(bact)^(1/4)

#Test for equal dispersion on fourth root transformed data 
frdist<-vegdist((froot), method="bray")
#Check that the treatment groups AREN'T significantly different to each other
TukeyHSD(betadisper (frdist,metadata$body.site))
#Check that the time groups AREN'T significantly different to each other
TukeyHSD(betadisper (frdist,metadata$year))
#Check that the subject groups AREN'T significantly different to each other
TukeyHSD(betadisper (frdist,metadata$subject))
```

### Do the PERMANOVA

### Do the PERMANOVA

We will run PERMANOVA on the transformed data, as it improved the dispersions. The formula contains a very important extra term: depth. This is a confounding variable that need to be taken into account.

We are using the updated function *adonis2*, with the argument *by = "margin"*. What this means is that the order of terms doesn't affect the *P* value shown for each one. This was not the case for the older function *adonis*. For more information on the topic, I recommend https://environmentalcomputing.net/statistics/linear-models/linear-regression/interpret-lm-coeffs/ and https://stats.stackexchange.com/questions/20452/how-to-interpret-type-i-type-ii-and-type-iii-anova-and-manova.

If an interaction is significant, *adonis2* will not show results for the main effects involved in that interaction.

One last note: every time you run PERMANOVA, expect to get slightly different *P* values. This is because there is a random element in how they are estimated.

```{r}
pma<-adonis2(froot ~ depth + Location + Sample.details + 
            Location:Sample.details,
            method="bray", data=metadata, by="margin", permutations=999)
#Look at results
pma
```

If you don't have a significant interaction, use the table above to get the P value for Sample.details, and run PERMANOVA again to get the P-value for site by listing site as the last of the main effects.

If you do have a significant interaction, there is no point in interpreting the main effects separately.

### Pairwise tests for PERMANOVA

Often you want to know not just whether a particular factor is significant, but where the differences occurs. For an ordinary ANOVA you would use something like a Tukey's post-hoc test, but for PERMANOVA we will have to code it manually. **DON'T PANIC!** The following code will walk you through the process with minimal pain.

#### Scenario 1: comparing the levels of a main effect.

The simplest type of pairwise comparison is comparing the levels of a main effect (a term not involved in an interaction). For example, you may want to know if treatment A differs from treatment B, and whether treatment B differs from treatment C.

Start by running the next section **without changing anything** to create the function we will use.

```{r}
main.effect.pairwise<-function(data, metadata, pair.var, other.terms=NULL, dist.method="bray", P.adj="BH", ...){
  modlist<-list()
  for(i in levels(metadata[,pair.var])[-1]){
    for(k in head(levels(metadata[,pair.var]), -1)){
      if (i != k){
        #Create a list containing subsets of the response variables in  pairwise combinations
        samp.pairs<-data[metadata[,pair.var] %in% c(i,k),]
        #Create a list containing subsets of the metadata 
        met.pairs<-metadata[metadata[,pair.var] %in% c(i,k),]
        #Do the comparisons
        if (!is.null(other.terms)){
          form<-as.formula(paste0("samp.pairs ~ ", pair.var, " + ", other.terms))
        } else
          form<-as.formula(paste0("samp.pairs ~ ", pair.var))
        pma.mod<-adonis2(form, data=met.pairs, method = dist.method, by="margin", ...)
        modlist[[paste(i,k,sep=".")]]<-pma.mod
        #Dispersion tests for  the pairwise combination datasets
        disp.model<-anova(betadisper(vegdist(samp.pairs, method=dist.method), met.pairs[,pair.var]))
        if (disp.model$`Pr(>F)`[!is.na(disp.model$`Pr(>F)`)] < 0.05){    print(paste("Warning: unequal dispersion for", i, "-", k, disp.model$`Pr(>F)`[!is.na(disp.model$`Pr(>F)`)], sep=" "))
        }
      }
    }
  }
  
  full.tab<-as.data.frame(t(sapply(modlist, function(x) x[pair.var,])))
  full.tab$P.adjusted <- p.adjust(full.tab$`Pr(>F)`, method = P.adj)
  
  return(full.tab)
}
```

Now you can actually run the comparisons! In this example, you're comparing pairs of treatments. Replace these with the column names of the variables you want to compare. The default distance method is "bray", but you can use any of the distances accepted by *adonis2*. The default method for adjusting *P* values for multiple comparisons is Benjamini-Hochberg, but you can use any method accepted by *p.adjust*. Any other arguments you wish to pass to *adonis2* can be added to the end.

If there are other terms you wish to include in the model (i.e. control for), you can add the argument *other.terms*. For example, if you wanted to control for the effect of temperature and genotype, you would add *other.terms="Temperature + Genotype"*.

```{r}
mfx.pairs<-main.effect.pairwise(data=bact, metadata=metadata, pair.var="Sample.details", dist.method="bray", P.adj="BH")
mfx.pairs
```

#### Scenario 2: an interaction where you are only interested in some of the possible comparisons.

Sometimes you want to explore an interaction only by comparing levels of one factor within levels of an interacting factor, rather than doing all possible comparisons. For example, you may want to know if treatment A differs from treatment B at time 2, but you're not interested in whether treatment A at time 1 differs from treatment B at time 2.

Start by running the next section **without changing anything** to create the function we will use.

```{r}
select.interaction.pairwise<-function(data, metadata, pair.var, single.var, other.terms=NULL, dist.method="bray", P.adj="BH", ...){
  modlist<-list()
  for(i in levels(metadata[,pair.var])){
    for(j in levels(metadata[,single.var])){
      for(k in levels(metadata[,pair.var])){
        if (grep(k, levels(metadata[,pair.var])) > grep(i, levels(metadata[,pair.var]))){
          print(paste(i,k, j))
          #Create a list containing subsets of the response variables in  pairwise combinations
          samp.pairs<-data[metadata[,pair.var] %in% c(i,k) & metadata[,single.var]==j,]
          #Create a list containing subsets of the metadata 
          met.pairs<-metadata[metadata[,pair.var] %in% c(i,k) & metadata[,single.var]==j,]
          #Do the comparisons
          if (!is.null(other.terms)){
            form<-as.formula(paste0("samp.pairs ~ ", pair.var, " + ", other.terms))
          } else
            form<-as.formula(paste0("samp.pairs ~ ", pair.var))
          pma.mod<-adonis2(form, data=met.pairs, method = dist.method, by="margin", ...)
          modlist[[paste(i,k,j,sep=".")]]<-pma.mod
          #Dispersion tests for  the pairwise combination datasets
          disp.model<-anova(betadisper(vegdist(samp.pairs, method=dist.method), met.pairs[,pair.var]))
          if (disp.model$`Pr(>F)`[!is.na(disp.model$`Pr(>F)`)] < 0.05){    print(paste("Warning: unequal dispersion for", i, "-", k, "-", j, disp.model$`Pr(>F)`[!is.na(disp.model$`Pr(>F)`)], sep=" "))
          }  
        }
      }
    }
  }
  
  full.tab<-as.data.frame(t(sapply(modlist, function(x) x[pair.var,])))
  full.tab$P.adjusted <- p.adjust(full.tab$`Pr(>F)`, method = P.adj)
  
  return(full.tab)
}
```

Now you can actually run the comparisons! In this example, you're comparing pairs of treatments within each time point. Replace these with the column names of the variables you want to compare. The default distance method is "bray", but you can use any of the distances accepted by *adonis*. The default method for adjusting *P* values for multiple comparisons is Benjamini-Hochberg, but you can use any method accepted by *p.adjust*. Any other arguments you wish to pass to *adonis* can be added to the end.

If there are other terms you wish to include in the model (i.e. control for), you can add the argument *other.terms*. For example, if you wanted to control for the effect of temperature and genotype, you would add *other.terms="Temperature + Genotype"*.

```{r}
select.pairs<-select.interaction.pairwise(data=bact, metadata=metadata, pair.var="Location", single.var="Sample.details", other.terms="depth", dist.method="bray", P.adj="BH")
select.pairs
```

#### Scenario 3: an interaction where you are interested in all the possible comparisons.

Start by running the next section **without changing anything** to create the function we will use.

```{r}
all.interaction.pairwise<-function(data, metadata, pair.var1, pair.var2, other.terms=NULL, dist.method="bray", P.adj="BH", ...){
  modlist<-list()
  metadata$interaction<-paste0(metadata[,pair.var1], metadata[,pair.var2])
for(i in levels(metadata[,pair.var1])){
  for(j in levels(metadata[,pair.var2])){
    for(k in levels(metadata[,pair.var1])){
        for(l in levels(metadata[,pair.var2])){
      if (!(i==k & j==l) & grep(k, levels(metadata[,pair.var1])) >= grep(i, levels(metadata[,pair.var1])) &
          grep(l, levels(metadata[,pair.var2])) >= grep(j, levels(metadata[,pair.var2]))) {
        print(paste(i,j,k,l))
    #Create a list containing subsets of the response variables in  pairwise combinations
    samp.pairs<-rbind(data[metadata[,pair.var1]==i & metadata[,pair.var2]==j,],
             data[metadata[,pair.var1]==k & metadata[,pair.var2]==l,])
    #Create a list containing subsets of the metadata 
    met.pairs<-rbind(metadata[metadata[,pair.var1]==i & metadata[,pair.var2]==j,],
             metadata[metadata[,pair.var1]==k & metadata[,pair.var2]==l,])
    #Do the comparisons
    if (!is.null(other.terms)){
    form<-as.formula(paste0("samp.pairs ~ interaction + ", other.terms))
    } else
      form<-as.formula(paste0("samp.pairs ~ interaction"))
    pma.mod<-adonis2(form, data=met.pairs, method = dist.method, by="margin", ...)
    modlist[[paste(i,k,j,l,sep=".")]]<-pma.mod
    #Dispersion tests for  the pairwise combination datasets
    disp.model<-anova(betadisper(vegdist(samp.pairs, method=dist.method), met.pairs[,"interaction"]))
  if (disp.model$`Pr(>F)`[!is.na(disp.model$`Pr(>F)`)] < 0.05){    print(paste("Warning: unequal dispersion for", i, "-", k, "-", j, "-", l, disp.model$`Pr(>F)`[!is.na(disp.model$`Pr(>F)`)], sep=" "))
   }  
   }
  }
 }
}
}
full.tab<-as.data.frame(t(sapply(modlist, function(x) x["interaction",])))
full.tab$P.adjusted <- p.adjust(full.tab$`Pr(>F)`, method = P.adj)

return(full.tab)
}

```

Now you can actually run the comparisons! In this example, you're comparing pairs of treatments within each time point. Replace these with the column names of the variables you want to compare. The default distance method is "bray", but you can use any of the distances accepted by *adonis*. The default method for adjusting *P* values for multiple comparisons is Benjamini-Hochberg, but you can use any method accepted by *p.adjust*. Any other arguments you wish to pass to *adonis* can be added to the end.

If there are other terms you wish to include in the model (i.e. control for), you can add the argument *other.terms*. For example, if you wanted to control for the effect of temperature and genotype, you would add *other.terms="Temperature + Genotype"*.

```{r}
all.pairs<-all.interaction.pairwise(data=bact, metadata=metadata, pair.var1="Location", pair.var2="Sample.details", other.terms = "depth", dist.method="bray", P.adj="BH")

#Look at the results
all.pairs
all.pairs[all.pairs$P.adjusted<0.05,]
```

# Metacoder plots for 16S rRNA gene data

Metacoder produces beautiful plots (heat trees) of taxonomic data. Documentation is available at https://grunwaldlab.github.io/metacoder_documentation/index.html

## Prepare the data for metacoder

```{r}
#Turn the data into proportions
#Turn the data around
bactprop<-t(bact)
#Take proportions and make into percentages
bactprop<-as.data.frame(prop.table(as.matrix(bactprop)))*100
```

## Read in the taxonomy information

```{r}
tax<-read_qza("taxonomy.qza")

#Add taxonomy column to the bactprop data frame
bactprop$taxonomy<-tax$data$Taxon[match(rownames(bactprop), tax$data$Feature.ID)]
#Tidy up names that haven't been formatted properly
bactprop$taxonomy<-gsub("\\.__", "", bactprop$taxonomy)
bactprop$taxonomy<-gsub("\\.([a-z]__)", ";\\1", bactprop$taxonomy)
```

## Create metacoder (taxmap) object

```{r}
#Remove reads unassigned at the Kingdom level
bactprop<-bactprop[grep("^Unassigned", bactprop$taxonomy, invert = T),]


#Create taxmap object
bacttaxmap <- parse_tax_data(bactprop,
                             class_cols = "taxonomy", # The column in the input table
                             class_sep = ";", # What each taxon is separated by
                             class_regex = "([a-z])__(.*)$",
                             class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
bacttaxmap

#Remove taxa without names
bacttaxmap <- filter_taxa(bacttaxmap, taxon_names != "")
#Remove Archaea
bacttaxmap <- filter_taxa(bacttaxmap, taxon_names == "Archaea", subtaxa = TRUE, invert=TRUE)

#'Collapse' the taxonomy to family level
bacttaxmap <- filter_taxa(bacttaxmap, taxon_ranks == "f", supertaxa = TRUE)
```

## Plot a 'first attempt' heat tree

```{r}
heat_tree(bacttaxmap,
          node_label = taxon_names,
          node_size = n_obs,
          node_color = n_obs)
```

## Plot heat trees for specific sample groups

Calculate read numbers for each taxon for each site.

```{r}
bacttaxmap$data$tax_abund <- calc_taxon_abund(bacttaxmap, "tax_data",
                                              cols = metadata$index,
                                              groups = metadata$body.site)
bacttaxmap$data$sum_abund <- rowSums(bacttaxmap$data$tax_abund[,-1])

#Plot a heat tree just for gut samples
bacttaxmap %>%
  filter_taxa(bacttaxmap$data$tax_abund$gut > 0) %>%
  heat_tree(
    node_label = taxon_names,
    node_size = gut,
    node_color = gut,
    initial_layout = "re", layout = "da",
    title = "Gut read depth",
    node_color_axis_label = "Sum of reads",
    node_size_axis_label = "Number of OTUs")
```

## Differential heat trees for body site

```{r metacoder-matrix}
#Calculate the mean difference between treatments for each compound, normalised to abundance
bacttaxmap$data$diff_table <- compare_groups(bacttaxmap,
                                             data = "tax_abund",
                                             cols = levels(metadata$body.site),
                                             groups = levels(metadata$body.site),
                                             func = function(abund1, abund2) 
                                               list(mean_diff = (mean(abund1) - mean(abund2)) / mean (abund1)))

#Put comparisons in the correct order
bacttaxmap$data$diff_table<-bacttaxmap$data$diff_table[with(bacttaxmap$data$diff_table, order(treatment_1,treatment_2)),]

#Draw heat tree matrix
bacttaxmap %>%
  heat_tree_matrix(data = "diff_table",
                   seed = 47,
                   node_size = sum_abund, #Use total abundance
                   node_label =  ifelse(sum_abund >= 1, as.character(taxon_names), NA),
                   node_label_size_trans = "log10 area",
                   node_label_size_range = c(0.02, 0.03),
                   node_color = mean_diff,
                   node_color_range = c("#FF2A00","#DDDDDD","#002AFF"),
                   node_color_trans = "area",
                   node_color_interval = c(-6, 6),
                   edge_color_interval = c(-6, 6),
                   row_label_color = "#002AFF",
                   col_label_color = "#DD2A00",
                   overlap_avoidance = 0.73,
                   node_size_axis_label = "Relative abundance",
                   node_color_axis_label = "Mean difference",
                   key_size = 0.78)
```

# RandomForests

Random forests are a machine learning technique used to predict whether data can be correctly classified into specific groups. (It is also possible to use random forests to predict values on a continuous scale, but for simplicity that isn't covered here.)
There is a clear and user-friendly  guide to random forests here:
https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm.
For (a lot) more technical detail, the original paper describing method is here: https://link.springer.com/article/10.1023%2FA%3A1010933404324.

#### Optimise the parameters

First, you need to optimise the settings for your data. Start by choosing the number of trees - more trees means more accuracy but also more computation time. Start by running a RF with a very large number of trees (100 000 in the code below). The resultant graph should show an exponentially decreasing error rate. Choose the approximate point on the x axis where the error has reached its asymptote, and then round up. This is the number of trees you will use.

```{r}
plot(randomForest(bact, y=metadata$body.site, proximity = F, keep.forest = F, ntree=100000))
```

Now you need to optimise the number of variables (*m*) used in building the trees. Substitute your choice of tree number for *ntree* in the code below. This time the error rate is likely to show a U-shape. Again, choose the value of x where the error rate is at its minimum. This is the value of *m* you will use.

```{r}
cv.treat<-rfcv(bact, metadata$body.site, cv.fold = 100, ntree = 5000)
plot(cv.treat$error.cv~cv.treat$n.var, type = "S", xlab = "No. variables", ylab = "Error rate", ylim = c(0,1), main = "RFCV")
```

## Run the random forest

Run the random forest (substitute your choice of tree number for ntree and your choice of *m* for mtry). The most important result from a RF is the confusion matrix. The rows represent the classes, and the columns represent the predicted classes. The matrix shows how many replicates from each true class were placed in each predicted class. The last column shows the proportion of replicates which were not correctly classified. 

```{r}
#Run a random forest classification for surface
siteRF<-randomForest(bact, y=metadata$body.site, importance = T, proximity = T, mtry = 25, ntree=50000)
#Look at how accurately the classification worked
siteRF$confusion

#Do an ordination plot based on the random forest
RFplot.surface<-MDSplot(siteRF, metadata$body.site, palette = cbfcol, xlim=c(-1,1), ylim=c(-1,1))
#Add 95% confidence intervals
ordiellipse(RFplot.surface, metadata$body.site,conf=0.95,col=cbfcol,lwd=2,draw="lines")
#Add a legend
legend("topleft",  c(levels(metadata$body.site)), pch=c(20), bty="n", cex=0.8, col=cbfcol)
```

## Find which variables were most important in the classification
The beauty of random forests is that they tell you not only how good the classification was, but also which variables (OTUs in the case of 16S data) were most important in determining the classification. This is calculated by randomly shuffling the values for each variable in turn, and comparing how much worse the error rate was for the shuffled data compared to the real data. The difference is the mean decrease in accuracy (MDA) for that variable. The random forest also calculates the Gini for each variable, which is a quick-and-dirty approximation of the MDA.
**N.B.** Note that importance scores are only comparable within, not between, random forests.

```{r}
#Plot two measures of which variables were most important
varImpPlot(siteRF, label=NA)
#If you want to, look at the whole list
importance(siteRF)
```

### Choose a subset of the most important variables

Often you will want to pick out the most important variables to explore them further. How many you choose is up to you, and there is no right answer. In looking at the graphs above, you can often see a 'break point' where there a group of variables with much higher MDAs than the rest. This example chooses the top 20 variables.

```{r}
#Extract the importance scores as a data frame
site.top.vars<-as.data.frame(importance(siteRF))
#Sort them in descending order of MDA and take the top 20
site.top.vars<-site.top.vars[order(site.top.vars$MeanDecreaseAccuracy, decreasing = T)[1:20],]

#See what they are!
tax$data$Taxon[match(rownames(site.top.vars), tax$data$Feature.ID)]
```
