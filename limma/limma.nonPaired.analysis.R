#/**
#* ##############################################################
#* Limma analysis
#* 2 group comparison
#* Created by [Varshna Goelela]
#* Updated: 2011 June 22
#* ##############################################################
#**/

## load source scripts
source("/Users/varshna/Documents/TNO/R scripts/limma/install.load.lib.R")
source("/Users/varshna/Documents/TNO/R scripts/limma/dateTime.R")

# variables and lists
#=====================
#choose a species:
sChoices <- c('Human','Mouse','Rat')
species <- select.list(sChoices)


# name study
ns  <- 'VetteKip'

#libs
lib.db <- paste('lumi',species,'All.db',sep='')
lib.mapping <- paste('lumi',species,'IDMapping',sep='')

man <- c( paste( 'lumi', species,'IDMapping', sep=''),
          paste( 'lumi', species,'All.db', sep=''),
          "limma",
          "ggplot2",
          "qvalue", "statmod"
          )

#* ###############################################
# Veranderen in fuctie
#* ###############################################

# load all the needed R library packages
loadPackages(man)

#folder containing Rdata files
indir <- "/Users/varshna/Documents/TNO/"

#select normData Rdata file object
Rdata <- c(list.files(indir, pattern = 'Rdata'))

#load normData R object
load(paste(indir, Rdata, sep=''))

#extract eset matrix from lumiBatch object file
data = exprs(normData)
#get list of nuIDs
nuIDs <- rownames(data)


#list text files in desc. folder
desc.table <- c(list.files(indir, pattern = 'description'))

#load desc file
descFile = paste(indir, desc.table[2], sep = "")
description <- read.table(descFile,
                          header=T,  
                          stringsAsFactors = F,
                          sep='\t',
                          quote="")



#check description file
#which column is the sample name
hdr = select.list(names(description), title = "Select the column with the sample names:")

#Match sampleNames from datafile with first column from description file
file_order <- match(description[,hdr],sampleNames(normData))
#Check on NA values in file_order
if(sum(is.na(file_order)) > 0) 
  stop("Assigned array names in raw data file and file names in description file do not match")
#Check if each value in file_order is unique
if(length(unique(file_order)) < length(file_order)) 
  stop("Assigned file names in description file are not unique")
#Check if length values description file is the same as unique length values description file
if(length(description[,hdr]) != length(unique(description[,hdr])) ) 
  stop("Assigned sampleNames are not unique")

cat("..::..::..\n", 
    "DISCRIPTION FILE OK!\n", sep="")

#reorder description 2
description = description[order(description[,hdr]),]



# specify the Subject type for each column. 
# collumns with the same name will be merged
Subjects <- factor(description$subjectNR)

# specify the Treatment type for each column. 
# collumns with the same name will be merged
Treatment <- factor(description$subGroupName)

#specify each timePoint for every sample
# collumns with the same name will be merged
groupTime <- factor(description$groupTime)

#determine levels groupTime
lev=levels(groupTime)



# Name comparison
contrast.fit <- "HFpio-HF"

if (require(limma)) {
  # create design matrix
  design <- model.matrix(~0+Treatment)
  colnames(design) <- lev
  rownames(design) <- description$sampleName
  fit <- lmFit(data, design) 
  #compare 2 groups
  cont.matrix <- makeContrasts(contrast.fit, levels=design)
  fit2 <- contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2)
  #Add gene symbols to gene properties
  if ( require(lib.db,  character.only = T) & require(annotate)) {
    
    symbol   <- unlist(lookUp(fit2$genes$ID, paste( 'lumi', species,'All.db', sep=''), "SYMBOL"))
    geneName <- unlist(lookUp(fit2$genes$ID, paste( 'lumi', species,'All.db', sep=''), "GENENAME")) 
    entrezID <- unlist(lookUp(fit2$genes$ID, paste( 'lumi', species,'All.db', sep=''), "ENTREZID"))
    accNum   <- unlist(lookUp(fit2$genes$ID, paste( 'lumi', species,'All.db', sep=''), "ACCNUM"))
    sequence <- unlist(id2seq(fit2$genes$ID))
    ilmnID   <- unlist(nuID2probeID(fit2$genes$ID, lib.mapping, species))
    
    #create dataframe with added meta data information
    fit2$genes   <- data.frame( nuID         = nuIDs,
                                probeID      = ilmnID,
                                geneSymbol   = symbol, 
                                geneName     = geneName, 
                                entrezID     = entrezID,
                                accessionNum = accNum,
                                Sequence     = sequence, stringsAsFactors = FALSE )                              
  }
  
  # Create a total results tabel with all the probes containg a pval below 1 and 
  # a log2 fold change above 0 and a infinite nr of probes.
  
  pval = 1      #cutoff value for adjusted p-values
  lfc  = 0      #cutoff value for log2-fold-change
  nr   = Inf    #maximum number of probes to list
  
  results = topTable(fit2, 
                     coef    = contrast.fit,
                     adjust  = 'fdr',
                     lfc     = lfc,
                     p.value = pval,
                     number  = nr)
}

#create a column with the Fold Change calculated from the logFC
results$FC <- 2^results$logF
#create a column with the calculated qvalues from the p.value
qobj <- qvalue(results$P.Value)
results$q.Value <- qobj$qvalues

#Write results in tab-delimeted table

results_fn = paste(sub("input", "output", indir),paste('LIMMA', ns, "nonPaired", contrast.fit, x=dateTime(), 'txt', sep='.'), sep="" )
write.table(results,  
            file      = results_fn,
            row.names = FALSE, 
            quote     = FALSE, 
            sep       = '\t')
