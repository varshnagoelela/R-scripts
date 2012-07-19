#/**
#* ##############################################################
#* Create cutt-off
#* Created by [Varshna Goelela]
#* Updated: 2011 July 11
#* ##############################################################
#**/
#empty all objects and memory
rm(list=ls(all=TRUE))

## load source scripts
source( "/Users/varshna/Documents/TNO/R scripts/install.load.lib.R")

#standard variabeles:
ns='RosiPio' #general name of the study

dateTIME= format(Sys.time(), "%Y%m%d_%H%M%S")


# create a list with all mandatory packages
man <- c('ggplot2', 'limma')

# load all the needed R library packages
loadPackages(man)

indir <- "/Users/varshna/Documents/TNO/"
setwd(indir)
getwd() #check if the directory is correct

#list files in dir
infiles = list.files(indir, pattern = 'txt')

#CREATE empty matrix containing nr of rows == nr of infiles
rowNMS = matrix ('NA', nrow= length(infiles) )
#rename colname of the empty matrix
colnames (rowNMS) = 'limmaComparisonFiles'

#rename rownames into the comparison in file names of infiles
for (i in 1:length(infiles)) {
  split_fn = unlist(strsplit(infiles[i], '\\.'))
  lbl = split_fn[3]
  rowNMS[i, 1] = lbl 
}


########################################################  

# create list with column names for datamatrix
colNMS <- c("FDR_0.01FC_1.5",
            "FDR_0.01",
            "FDR_0.05FC_1.5",
            "FDR_0.05",
            "FDR_0.1",
            "Pval_0.01FC_2",  
            "Pval_0.01FC_1.5",	
            "Pval_0.01",
            "Pval_0.05FC_2",	
            "Pval_0.05FC_1.5",	
            "Pval_0.05" )

# create a empty datamatrix containing 10 columns and 
# as many rows as there are files in infiles
cutoffM <- matrix(0, length(infiles), 11, dimnames =list(infiles,colNMS))
#rownames(cutoffM) <- infiles

# for each file in list infiles
for (fn in infiles) {
  
  # read table fn
  sample = read.file(fn)
  
  # save how many probeids meet the cutoff criteria
  #CUTOFF_1: FDR_0.01FC_1.5:
  c01= which ( ( (sample$adj.P.Val < 0.01) & (sample$logFC > 0.5849625)  ) | 
               ( (sample$adj.P.Val < 0.01) & (sample$logFC < -0.5849625) ) )
  cutoffM[fn, 1] = length(c01)

  #CUTOFF_2: FDR_0.01:              
  c02 = which (sample$adj.P.Val < 0.01)
  cutoffM[fn, 2] = length(c02)
                
  #CUTOFF_3: FDR_0.05FC_1.5:  
  c03= which ( ( (sample$adj.P.Val < 0.05) & (sample$logFC > 0.5849625)  ) | 
               ( (sample$adj.P.Val < 0.05) & (sample$logFC < -0.5849625) ) )
  cutoffM[fn, 3] = length(c03)
  
  #CUTOFF_4: FDR_0.05:
  c04 = which (sample$adj.P.Val < 0.05)
  cutoffM[fn, 4] = length(c04)
  
  #CUTOFF_5: FDR_0.1:
  c05 = which (sample$adj.P.Val < 0.1)
  cutoffM[fn, 5] = length(c05)

  #CUTOFF_6: Pval_0.01FC_2:
  c06= which ( ( (sample$P.Value < 0.01) & (sample$logFC > 1)  ) | 
               ( (sample$P.Value < 0.01) & (sample$logFC < -1) ) )
  cutoffM[fn, 6] = length(c06)
  
  #CUTOFF_7: Pval_0.01FC_1.5:
  c07= which ( ( (sample$P.Value < 0.01) & (sample$logFC > 0.5849625)  ) | 
               ( (sample$P.Value < 0.01) & (sample$logFC < -0.5849625) ) )
  cutoffM[fn, 7] = length(c07)
  
  #CUTOFF_8: Pval_0.01:
  c08 = which (sample$P.Value < 0.01)
  cutoffM[fn, 8] = length(c08)
  
  #CUTOFF_9: Pval_0.05FC_2:
  c09= which ( ( (sample$P.Value < 0.05) & (sample$logFC > 1)  ) | 
               ( (sample$P.Value < 0.05) & (sample$logFC < -1) ) )
  cutoffM[fn, 9] = length(c09)
                        
  #CUTOFF_10: Pval_0.05FC_1.5: 
  c10= which ( ( (sample$P.Value < 0.05) & (sample$logFC > 0.5849625)  ) | 
               ( (sample$P.Value < 0.05) & (sample$logFC < -0.5849625) ) )
  cutoffM[fn, 10] = length(c10)
                                              
  #CUTOFF_11: Pval_0.05:
  c11 = which (sample$P.Value < 0.05)
  cutoffM[fn, 11] = length(c11)
}
                      
cutoffM = cbind(rowNMS, cutoffM)

new.indir <- paste(indir, 'Cutoff', dateTIME,  sep='/')
dir.create(new.indir)
setwd(new.indir)

# create a file name and save the data.frame as a tab-delimited text file
file.nm = paste(ns, 'cutoff', dateTIME,'txt', sep='.')
write.table(cutoffM,
            file = file.nm,
            quote= FALSE,
            sep='\t',
            row.names= F,
            col.names= T )

file.nm = paste(ns, 'cutoff', dateTIME,'xls', sep='.')
write.table(cutoffM,
            file = file.nm,
            quote= FALSE,
            sep='\t',
            row.names= F,
            col.names= T )

# redirect back to original workspace
setwd(indir)
        