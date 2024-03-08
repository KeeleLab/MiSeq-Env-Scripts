# 2024_03_07
# v2v3a4b7 analysis code - If you are reading this, please note that Charles Goodman will likely be updating this 
# to be performed using read mapping rather than positional searching in R. 

# load runInfo form used for analysis
# You will need two run info forms, one for V2V3 primers and one for A4B7 primers
# If you only have one form with V2 and A4B7 primers in it, please split it into two forms before proceeding

# please note.. both run info forms are meant to have the same samples in each run info form
# if you are missing a sample in either run, the final formatting will be incorrect and you will need to manually
# generate the output file
# Ideally samples will contain the same P5 primer, although this doesn't necessarily have to be
# the case as long as the samples remain in the same order

# read in common functions
library(openxlsx)
library(dplyr)
read_all_sheets = function(xlsxFile) {
  tryCatch(
    for (i in xlsxFile){
      sheet_names = openxlsx::getSheetNames(xlsxFile)                               # Get sheet names
      sheet_list = as.list(rep(NA, length(sheet_names)))                            # Initialize an empty list with space for all sheets
      names(sheet_list) = sheet_names                                               # Name list items with sheet names
      # For every sheet in the list of names:
      for (sn in sheet_names) {
        # add the actual data to the list of sheets
        sheet_list[[sn]] = openxlsx::read.xlsx(xlsxFile, sheet=sn, rowNames = F, colNames = T, skipEmptyRows = FALSE, skipEmptyCols = FALSE)
      }
      return(sheet_list)
    }, error = function(c){
      c$message <- paste0(c$message, "\n \nThe script encountered a problem opening the following file\n\n'", i, "' \n\nPlease try opening the file in Excel to ensure that it works. If the file has the ~ character in front of its name, the file is currently open.\n
If this is the case, try closing the file and rerunning the script. \n\n")
      stop(c)
    })
}
## Provide Run Info form paths
runInfoFormV2V3 <- "pathToRunInfoForm1.xlsx"
runInfoFormA4B7 <- "pathToRunInfoForm2.xlsx"
runInfoList <- list(runInfoFormV2V3, runInfoFormA4B7)
anyFalse <- sapply(runInfoList, is.null)
if(sum(anyFalse) > 0){
  runInfoList <- runInfoList[-which(anyFalse == T)]
  stop("One of the run info forms does not exist. Please double check your forms before continuing. Forms should have identical sample names, and primer numbers associated with them. 
       For example, sample 1 should have primer V2V3.P5.1 and A4B7.P5.1. Otherwise the final formatting will be incorrect.")
}
runInfoNames <- list()
runInfoList <- lapply(runInfoList, function(file){
  runInfo <- read_all_sheets(file)[[1]]
  names1 <- runInfo$Animal[!is.na(runInfo$Animal)]
  names2 <- runInfo$Sample[!is.na(runInfo$Sample)]
  if(is.character(runInfo$Date)){
    runInfo$Date <- as.character(as.Date(runInfo$Date, tryFormats = c("%m/%d/%y", "%Y-%m-%d")))
  }else if(is.numeric(runInfo$Date)){
    runInfo$Date <- as.character(as.Date(runInfo$Date, origin = "1899-12-30"))
  }
  names3 <- runInfo$Date[!is.na(runInfo$Date)]
  names4 <- gsub(".*([0-9]){1,2}$", "\\1", runInfo$Barcodes[!is.na(runInfo$Barcodes)])
  names <- paste(names1,names2,names3,names4)
  runInfoNames[[file]] <<- names
  return(runInfo)
})
## At this point, please check that runInfoNames[1] pairs with runInfoNames[2]. If it does not, DO NOT CONTINUE PAST
## THE FINAL FORMATTING STAGE. Manual adjustment will be required.
as.data.frame(cbind(runInfoNames[[1]], runInfoNames[[2]]))

if(length(runInfoNames) == 2){
  if(length(runInfoNames[[1]]) != length(runInfoNames[[2]])){
    stop("Run Info Forms Do Not Contain the Same Number of Samples - Manual Output Required. Manually need to combine stats values with Run Info Forms and assign names to concatDFs")
  }
  names <- runInfoNames[[1]]
}

runInfoList <- bind_rows(runInfoList)
names(runInfoList) <- c("Run.Number",
                        "Animal",
                        "Sample",
                        "Date",
                        "Barcodes",
                        "Input.TOTAL.PER.BARCODE",
                        "Input.TOTAL.PER.WELL",
                        "",
                        "(F.Barcode)",
                        "(cDNA)")
runInfoList$Barcodes <- toupper(runInfoList$Barcodes)
runInfoList <- runInfoList[,-1]
# remove any NA rows so that we can bind the run info forms
if(sum(is.na(runInfoList$Barcodes)) > 1){
  runInfoList <- runInfoList[-which(is.na(runInfoList$Barcodes)),]
}


# Need to primer match then find the v2 loop on R2

# 5' indexing tag used during sequencing

## For now we manually specify the primers to be used.. but could extract from runinfo form

# Primers used in V2V3 Run
V2V3primers <- c(1:4)
# Primers used in A4B7 Run 
a4b7primers <- c(9:12)


barcodesused <- V2V3primers

# Path to R1
reversefilename <- "pathToV2V3R1.fastq.gz"

# Path to R2y
forwardfilename <- "pathToV2V3R2.fastq.gz"

library(ShortRead)
library(stringr)
library(Biostrings)
library(tibble)
library(openxlsx)
library(dplyr)
# install.packages("qpcR")
library(qpcR)
library(stats)
## bug on some machines with "rgl" package loaded by qpcR, this circumvents that
options(rgl.useNULL= TRUE)
revfile <- readFastq(reversefilename)
# # Read in the reverse sequences
rev.seq <- as.character(sread(revfile))
# # And remove the rest of the file for memory reasons
rm(revfile)
gc()

# Trim adapters from R1s
adapter <- "NNNN"
rev.seq <- str_sub(rev.seq, nchar(adapter) + 1, nchar(rev.seq))

primer <- paste("P5.",1:40,sep="") #writes primer names, doesn't really matter if its VPX, INT, PBS, etc. 
barcodeseqs <- c("TAAGGCGA","CGTACTAG","AGGCAGAA","TCCTGAGC","GGACTCCT","TAGGCATG","CTCTCTAC","CAGAGAGG","GCTACGCT","CGAGGCTG",
                 "AAGAGGCA","GTAGAGGA","TAGATCGC","CTCTCTAT","TATCCTCT","AGAGTAGA","GTAAGGAG","ACTGCATA","AAGGAGTA","CTAAGCCT",
                 "TATAGCCT","ATAGAGGC","CCTATCCT","GGCTCTGA","AGGCGAAG","TAATCTTA","CAGGACGT","GTACTGAC","ATTACTCG",'TCCGGAGA',
                 "CGCTCATT","GAGATTCC","ATTCAGAA","GAATTCGT","CTGAAGCT","TAATGCGC","CGGCTATG","TCCGCGAA",'TCTCGCGC','AGCGATAG') #40 multiplexed barcodes

barcode.list <- cbind(primer, barcodeseqs) #this holds the primer number, sequence, inputs, and whether to keep 
# Pull the primer sequences from the R1s
barcodes <- str_sub(rev.seq, 1, 8)
barcodematches <- match(barcodes, barcode.list[,2])

# Table to hold summary statistics
stats <- data.frame(Primer = character(), "# of 5' Indexing Sequences" = numeric(),
                    "Reads matching V2 reference" = numeric(), stringsAsFactors = FALSE)

rm(rev.seq)
gc()

fwdfile <- readFastq(forwardfilename)
# Read in the forward sequences
fwd.seq <- as.character(sread(fwdfile))
# And remove the rest of the file because memory
rm(fwdfile)
gc()



#i = 6
V2translations <- list()
V2barcodes <- list()
for(i in barcodesused){
  print(paste("Primer", i, barcode.list[i,2]))
  
  # Find the reads matching to this primer
  splitmatch=(barcodematches==i)&(!is.na(barcodematches))
  
  # pull the length of the read after the primer
  
  # We only actually care about the forward reads for this one
  # Pull all the R2s for the reads that matched this primer
  fwd.seq.split.save <- fwd.seq[splitmatch]
  
  # Let's look for the reference sequence from Brandon
  
  ref <- "AACAAGAGCAAATGATAAGC"
  pos <- unlist(aregexec(ref, fwd.seq.split.save, max.distance = 5))
  # pos <- pos[pos > 0 & pos <= 150]
  print(paste("Matches:", length(pos[pos > 0])))
  
  # We want the 132 bases downstream of this.
  v2 <- substr(fwd.seq.split.save[pos > 0], pos[pos > 0] + 20, pos[pos > 0] + 20 + 132 - 1)
  v2 <- v2[which(nchar(v2) == 132)]
  # Make a table of the most common sequences
  table <- as.data.frame(table(v2))
  table$AA <- as.character(translate(DNAStringSet(table$v2), if.fuzzy.codon = "solve"))
  
  # And also a table of just the translations
  trans <- table[,2:3]
  trans <- aggregate(Freq ~ AA, trans, sum)
  
  # Discard singletons
  table <- table[(table$Freq > 1),]
  trans <- trans[(trans$Freq > 1),]
    
  # Calculate proportions
  table$proportion <- table$Freq/sum(table$Freq)
  trans$proportion <- trans$Freq/sum(trans$Freq)
  
  # Order counts by frequency
  table <- table[order(table$Freq, decreasing = TRUE),]
  trans <- trans[order(trans$Freq, decreasing = TRUE),]
  
  # Clean up tables a bit
  if(length(table[[1]]) > 0){
    table <- table[,c(1,3,2,4)]
    rownames(table) <- 1:length(table[[1]])
  }
  if(length(trans[[1]]) > 0){
    rownames(trans) <- 1:length(trans[[1]])
  }
  
  # Update our summary table
  stats[nrow(stats)+1,] <- c(paste0("V2V3.", barcode.list[i]), sum(splitmatch), length(pos[pos > 0]))

  # Need to make them the same length so we can output them in one file
  holdlen <- max(length(table[[1]]), length(trans[[1]]))
  
  # Actually making them one longer than the longest row so that nothing gets written over
  # (Fixing bug with Sirish method where the last row gets removed if the table is already the desired length)
  
  table[nrow(table) + (holdlen - nrow(table)) + 1,] <- NA
  trans[nrow(trans) + (holdlen - nrow(trans)) + 1,] <- NA
  
  # And now let's combine these into one dataframe
  table[[(ncol(table)) + 1]] <- " "
  hold <- data.frame(table, trans, stringsAsFactors = FALSE)
  
  colnames(hold) <- c("V2 Sequence", "V2 Translation", "V2 Count", "V2 Proportion", "", "V2 Translation", "V2 Count", "V2 Proportion")

  names(trans) <- c("V2 Sequence", "V2 Count", "V2 Proportion")
  V2translations[[paste0(i)]] <- trans[complete.cases(trans),]
  V2barcodes[[paste0(i)]] <- table
}
names(V2translations) <- barcodesused
names(V2barcodes)<-barcodesused

colnames(stats) <- c("Primer", "# of 5' Indexing Sequences", "# of Reads Matching V2 Reference")
V2Stats <- stats
write.csv(stats, "pathToOutputStatsFileForV2.csv", na = "", row.names = FALSE)

rm(fwd.seq)
rm(barcodes)
rm(barcodematches)
rm(fwd.seq.split.save)
rm(splitmatch)
gc()
###A4B7 setup ###
# Path to A4B7 R1
reversefilename <- "pathToA4B7R1.fastq.gz"

# Path to A4B7 R2
forwardfilename <- "pathToA4B7R2.fastq.gz"

revfile <- readFastq(reversefilename)
# Read in the reverse sequences
rev.seq <- as.character(sread(revfile))
rev.qual <- as.character(quality(revfile)@quality)
# And remove the rest of the file for memory reasons
rm(revfile)
gc()

fwdfile <- readFastq(forwardfilename)
# Read in the forward sequences
fwd.seq <- as.character(sread(fwdfile))
# And remove the rest of the file because memory
rm(fwdfile)
gc()

# Trim adapters from R1s
adapter <- "NNNN"
rev.seq <- str_sub(rev.seq, nchar(adapter) + 1, nchar(rev.seq))
gc()

primer <- paste("P5.",1:40,sep="") #writes primer names, doesn't really matter if its VPX, INT, PBS, etc. 
barcodeseqs <- c("TAAGGCGA","CGTACTAG","AGGCAGAA","TCCTGAGC","GGACTCCT","TAGGCATG","CTCTCTAC","CAGAGAGG","GCTACGCT","CGAGGCTG",
                 "AAGAGGCA","GTAGAGGA","TAGATCGC","CTCTCTAT","TATCCTCT","AGAGTAGA","GTAAGGAG","ACTGCATA","AAGGAGTA","CTAAGCCT",
                 "TATAGCCT","ATAGAGGC","CCTATCCT","GGCTCTGA","AGGCGAAG","TAATCTTA","CAGGACGT","GTACTGAC","ATTACTCG",'TCCGGAGA',
                 "CGCTCATT","GAGATTCC","ATTCAGAA","GAATTCGT","CTGAAGCT","TAATGCGC","CGGCTATG","TCCGCGAA",'TCTCGCGC','AGCGATAG') #40 multiplexed barcodes

barcode.list <- cbind(primer, barcodeseqs) #this holds the primer number, sequence, inputs, and whether to keep 
# Pull the primer sequences from the R1s
barcodes <- str_sub(rev.seq, 1, 8)
barcodematches <- match(barcodes, barcode.list[,2])

#### a4b7 analysis part
# Table to hold summary statistics
stats <- data.frame(Primer = character(), "# of 5' Indexing Sequences" = numeric(),
                    "MutA reads" = numeric(), "MutB reads" = numeric(),
                    "MutC reads" = numeric(), "MutD reads" = numeric(),
                    "MutE reads" = numeric(), stringsAsFactors = FALSE)

# Annotations for the first mutation
# mutAlist <- read.csv("file:///C:/Users/thorpeal/Documents/Keele Lab/121118 Mario Miseq prep/Mutation A annotations.csv", 
#                      stringsAsFactors = FALSE)

barcodesused <- a4b7primers

#i = 1
a4b7trans <- list()
a4b7barcodes <- list()
for(i in barcodesused){
  print(paste("Primer", i, barcode.list[i,2]))
  
  # Find the reads matching to this primer
  splitmatch=(barcodematches==i)&(!is.na(barcodematches))
  
  print(paste("Matches: ", sum(splitmatch)))
  
  # pull the length of the read after the primer
  rev.seq.split.save <- str_sub(rev.seq[splitmatch], nchar(barcode.list[i,2]) + 1, nchar(rev.seq[splitmatch]))
  
  # Pull all the R2s for the reads that matched this primer
  fwd.seq.split.save <- fwd.seq[splitmatch]
  
  # Let's look for the first mutation first...
  refA <- "AAAAGGAATTACGTGCCATGTCATATTAGACAAATA"
  # Want the reads reverse complemented
  rev.seq.fwd <- as.character(reverseComplement(DNAStringSet(rev.seq.split.save)))
  
  # Find the reference sequence
  posA <- unlist(aregexec(refA, rev.seq.fwd, max.distance = 5))

  print(paste("MutA:", length(posA[posA > 0])))
  # Now let's grab the 37 bases next to it...
  mutA <- substr(rev.seq.fwd[posA > 0], posA[posA > 0] - 37, posA[posA > 0] - 1)
  # And see what we find there
  tableA <- as.data.frame(table(mutA))
  tableA$mutA <- as.character(tableA$mutA)
  
  # Translate the sequence to amino acids to make it easier to interpret
  protein <- substr(tableA$mutA, 2, nchar(tableA$mutA[1]))
  protein <- as.character(translate(DNAStringSet(protein), if.fuzzy.codon = "solve"))
  tableA$AA <- protein
  
  # Match sequences to our annotations
  colnames(tableA)[1] <- "DNA.Sequence"
  # tableA <- merge(tableA, mutAlist, by = "DNA.Sequence", all.x = TRUE)
  
  # Let's do some rearranging to clean up this table and put the annotations on the side
  if(length(tableA[[1]]) > 0){
    # Order counts by frequency
    tableA <- tableA[order(tableA$Freq, decreasing = TRUE),]
    # Remove singletons
    tableA <- tableA[(tableA$Freq > 1),]
    # And calculate proportions
    tableA$proportion <- tableA$Freq/sum(tableA$Freq)
    # This is how many we found in the end
    sum(tableA$Freq)
    tableA <- tableA[,c(1,3,2,4)]
    # And give the column headers a bit more meaningful names
    colnames(tableA) <- c("MutA Sequence", "MutA Translation","Counts",  "Proportions")
  }else{
    tableA$proportions <- numeric(0)
  }
  
  # Next mutation, this one on the forward strand so I haven't found it with the app to compare
  refB <- "CTATAATTAGTTTAAATAAGTA"
  
  posB <- unlist(aregexec(refB, fwd.seq.split.save, max.distance = 3))
  print(paste("MutB:", length(posB[posB > 0])))
  # For this we only want the base directly before the reference, but let's pull the whole codon in case
  mutB <- substr(fwd.seq.split.save[posB > 0], posB[posB > 0] - 16, posB[posB > 0] + 16)
  tableB <- as.data.frame(table(mutB))
  tableB$AA <- as.character(translate(DNAStringSet(tableB$mutB), if.fuzzy.codon = "solve"))
  
  # Let's clean this up a bit
  if(length(tableB[[1]]) > 0 ){
    # I think just having the translation is fine for these, don't really need annotations
    # Discard singletons
    tableB <- tableB[(tableB$Freq > 1),]
    # Calculate proportions
    tableB$proportion <- tableB$Freq/sum(tableB$Freq)
    # Order counts by frequency
    tableB <- tableB[order(tableB$Freq, decreasing = TRUE),]
    tableB <- tableB[,c(1,3,2,4)]
    colnames(tableB) <- c("mutB Codon", "mutB Translation", "Counts", "Proportions")
  }else{
    tableB$Proportions <- numeric(0)
    tableB$`mutB Translation` <- numeric(0)
  }
  
  # Next mutation, also on forward strand. 25 nt reference sequence upstream
  
  refC <- "ACCATTGTCAAACATCCCAGGTATA"
  
  posC <- unlist(aregexec(refC, fwd.seq.split.save, max.distance = 3))
  print(paste("MutC:", length(posC[posC > 0])))
  
  # Swapped this around so now want the base just downstream of reference, taking the whole codon too
  mutC <- substr(fwd.seq.split.save[posC > 0], posC[posC > 0] + 25 - 7, posC[posC > 0] + 25 + 25)
  tableC <- as.data.frame(table(mutC))
  
  # Translate the codons
  tableC$AA <- as.character(translate(DNAStringSet(tableC$mutC), if.fuzzy.codon = "solve"))
  
  if(length(tableC[[1]]) > 0){
    # Remove singletons
    tableC <- tableC[(tableC$Freq > 1),]
    # Calculate proportions
    tableC$proportions <- tableC$Freq/sum(tableC$Freq)
    # Order counts by frequency
    tableC <- tableC[order(tableC$Freq, decreasing = TRUE),]
    # Clean up a bit
    tableC <- tableC[,c(1,3,2,4)]
    colnames(tableC) <- c("mutC Codon", "mutC Translation", "Counts", "Proportions")
  }else{
    tableC$Proportions <- numeric(0)
    tableC$'muB Trnaslation' <- numeric(0)
  }
  
  # Fourth region which is actually pretty close to the third (and also pretty close to the end of the fwd read)
  
  refD <- "ACTGATAAAATCAATTTGACGGCTCCT"
  
  posD <- unlist(aregexec(refD, fwd.seq.split.save, max.distance = 3))
  print(paste("MutD:", length(posD[posD > 0])))
  
  # We just want the one base at the end of the reference sequence, grabbing the whole codon though
  mutD <- substr(fwd.seq.split.save[posD > 0], posD[posD > 0] + 12, posD[posD > 0] + 27 + 2)
  tableD <- as.data.frame(table(mutD))
  # Translate the codons
  tableD$AA <- as.character(translate(DNAStringSet(tableD$mutD), if.fuzzy.codon = "solve"))
  
  if(length(tableD[[1]]) > 0){
    # Remove singletons
    tableD <- tableD[(tableD$Freq > 1),]
    # Calculate proportions
    tableD$proportion <- tableD$Freq/sum(tableD$Freq)
    # Order by frequency
    tableD <- tableD[order(tableD$Freq, decreasing = TRUE),]
    # Clean up a bit
    tableD <- tableD[,c(1,3,2,4)]
    colnames(tableD) <- c("mutD Codon", "MutD Translation", "Counts", "Proportions")
  }else{
    tableD$proportions <- numeric(0)
    tableD$translation <- numeric(0)
  }
  
  # Last mutation is on the reverse read, right after the primer
  refE <- "ACCAGTCTCATAGCAAACATAGATTGG"
  
  posE <- unlist(aregexec(refE, rev.seq.fwd, max.distance = 3))
  print(paste("MutE:", length(posE[posE > 0])))
  
  # We want the 25 nt downstream of this reference
  # Actually going to grab more to get the full codon
  mutE <- substr(rev.seq.fwd[posE > 0], posE[posE > 0] + 27 +3, posE[posE > 0] + 27 + 25 + 4)
  # Site 479 N->S mutation
  # mutE1 <- substr(rev.seq.fwd[posE > 0], posE[posE > 0] + 27 + 18, posE[posE > 0] + 27 + 18 + 2)
  # # Site 481 T->A mutation
  # mutE2 <- substr(rev.seq.fwd[posE > 0], posE[posE > 0] + 27 + 24, posE[posE > 0] + 27 + 24 + 2)
  
  tableE <- as.data.frame(table(mutE))
  # tableE1 <- as.data.frame(table(mutE1))
  # tableE2 <- as.data.frame(table(mutE2))
  # Translate to protein
  # Translate the codons
  tableE$AA <- as.character(translate(DNAStringSet(tableE$mutE), if.fuzzy.codon = "solve"))
  
  if(length(tableE[[1]]) > 0){
    # Remove singletons
    tableE <- tableE[(tableE$Freq > 1),]
    # Calculate proportions
    tableE$proportion <- tableE$Freq/sum(tableE$Freq)
    # Order by frequency
    tableE <- tableE[order(tableE$Freq, decreasing = TRUE),]
    # Clean up a bit
    tableE <- tableE[,c(1,3,2,4)]
    colnames(tableE) <- c("mutD Codon", "MutE Translation", "Counts", "Proportions")
  }else{
    tableE$proportions <- numeric(0)
    tableE$translation <- numeric(0)
  }
  
  stats[nrow(stats)+1,] <- c(paste0("A4B7.", barcode.list[i]), sum(splitmatch), length(posA[posA > 0]), length(posB[posB > 0]), 
                             length(posC[posC > 0]), length(posD[posD > 0]), length(posE[posE > 0]))
  
  # Adding a separate output file that's grouped by translation instead of nucleotide sequence
  transA <- tableA[,2:3]
  transA <- transA[complete.cases(transA),]
  if(nrow(transA) > 0){
    transA <- aggregate(Counts ~ `MutA Translation`, transA, sum)
    transA$Proportions <- transA$Counts/sum(transA$Counts)
    transA <- transA[order(transA$Counts, decreasing = TRUE),]
    rownames(transA) <- 1:length(transA[[1]])
  }else{
    transA$proportions <- numeric(0)
  }
  
  transB <- tableB[,2:3]
  transB <- transB[complete.cases(transB),]
  if(nrow(transB) > 0){
    transB <- aggregate(Counts ~ `mutB Translation`, transB, sum)
    transB$Proportions <- transB$Counts/sum(transB$Counts)
    transB <- transB[order(transB$Counts, decreasing = TRUE),]
    rownames(transB) <- 1:length(transB[[1]])
  }else{
    transB$proportions <- numeric(0)
    transB$Counts <- numeric(0)
    
  }
  
  transC <- tableC[,2:3]
  transC <- transC[complete.cases(transC),]
  if(nrow(transC) > 0){
    transC <- aggregate(Counts ~ `mutC Translation`, transC, sum)
    transC$Proportions <- transC$Counts/sum(transC$Counts)
    transC <- transC[order(transC$Counts, decreasing = TRUE),]
    rownames(transC) <- 1:length(transC[[1]])
  }else{
    transC$proportions <- numeric(0)
    transC$Counts <- numeric(0)
  }
  
  transD <- tableD[,2:3]
  transD <- transD[complete.cases(transD),]
  if(nrow(transD) > 0){
    transD <- aggregate(Counts ~ `MutD Translation`, transD, sum)
    transD$Proportions <- transD$Counts/sum(transD$Counts)
    transD <- transD[order(transD$Counts, decreasing = TRUE),]
    rownames(transD) <- 1:length(transD[[1]])
    transD <- transD[complete.cases(transD),]
  }else{
    transD$proportions <- numeric(0)
    transD$Counts <- numeric(0)
  }
  
  transE <- tableE[,2:3]
  transE <- transE[complete.cases(transE),]
  if(nrow(transE) > 0){
    transE <- aggregate(Counts ~ `MutE Translation`, transE, sum)
    transE$Proportions <- transE$Counts/sum(transE$Counts)
    transE <- transE[order(transE$Counts, decreasing = TRUE),]
    rownames(transE) <- 1:length(transE[[1]])
    transE <- transE[complete.cases(transE),]
  }else{
    transE$proportions <- numeric(0)
    transE$Counts <- numeric(0)
  }
  
  # So first of all we need to format the five mutation tables and make sure they're all the same length so we can combine them
  if (nrow(tableA) > 0){
    rownames(tableA) <- 1:length(tableA[[1]])
  }
  if(nrow(tableB) > 0){
    rownames(tableB) <- 1:length(tableB[[1]]) 
  }
  if(nrow(tableC) > 0){
    rownames(tableC) <- 1:length(tableC[[1]]) 
  }
  if(nrow(tableD) > 0){
    rownames(tableD) <- 1:length(tableD[[1]]) 
  }
  if(nrow(tableE) > 0){
    rownames(tableE) <- 1:length(tableE[[1]])
  }
  
  # Okay now that these are all formatted nicely, let's make them all the same length
  holdlen <- max(length(tableA[[1]]), length(tableB[[1]]), length(tableC[[1]]), length(tableD[[1]]),
              length(tableE[[1]]))
                # length(tableE1[[1]]), length(tableE2[[1]]))
  
  # Actually making them all one longer than the longest row so that nothing gets written over
  # (Fixing bug with Sirish method where the last row gets removed if the table is already the desired length)
  
  tableA[nrow(tableA) + (holdlen - nrow(tableA)) + 1,] <- NA
  tableB[nrow(tableB) + (holdlen - nrow(tableB)) + 1,] <- NA
  tableC[nrow(tableC) + (holdlen - nrow(tableC)) + 1,] <- NA
  tableD[nrow(tableD) + (holdlen - nrow(tableD)) + 1,] <- NA
  tableE[nrow(tableE) + (holdlen - nrow(tableE)) + 1,] <- NA
  # tableE1[nrow(tableE1) + (holdlen - nrow(tableE1)) + 1,] <- NA
  # tableE2[nrow(tableE2) + (holdlen - nrow(tableE2)) + 1,] <- NA
  
  # And now let's combine these all into one dataframe, agh...
  tableA[[(ncol(tableA)) + 1]] <- " "
  hold <- data.frame(tableA, tableB, stringsAsFactors = FALSE)
  hold[[ncol(hold) + 1]] <- NA
  
  hold <- data.frame(hold, tableC, stringsAsFactors = FALSE)
  hold[[ncol(hold) + 1]] <- NA
  
  hold <- data.frame(hold, tableD, stringsAsFactors = FALSE)
  hold[[ncol(hold) + 1]] <- NA
  
  hold <- data.frame(hold, tableE, stringsAsFactors = FALSE)
  hold[[ncol(hold) + 1]] <- NA
  
  colnames(hold) <- c("MutA Sequence", "Translation", "Count", "Proportion",
                      "", "MutB (297 T->A)", "Translation", "Count", "Proportion", "", "MutC (368 T->N)", "Translation", "Count", "Proportion",
                      "", "MutD (382 G->R)", "Translation", "Count", "Proportion", "", "MutE1 (479 N->S)", "Translation",
                      "Count", "Proportion", "" )
  
  
  # And now combine the translation tables into a separate file
  holdlen <- max(length(transA[[1]]), length(transB[[1]]), length(transC[[1]]), length(transD[[1]]), 
  length(transE[[1]]))
  
  transA[nrow(transA) + (holdlen - nrow(transA)) + 1,] <- NA
  transB[nrow(transB) + (holdlen - nrow(transB)) + 1,] <- NA
  transC[nrow(transC) + (holdlen - nrow(transC)) + 1,] <- NA
  transD[nrow(transD) + (holdlen - nrow(transD)) + 1,] <- NA
  transE[nrow(transE) + (holdlen - nrow(transE)) + 1,] <- NA
  
  transA[[(ncol(transA)) + 1]] <- " "
  trans <- data.frame(transA, transB, stringsAsFactors = FALSE)
  trans[[ncol(trans) + 1]] <- NA
  
  trans <- data.frame(trans, transC, stringsAsFactors = FALSE)
  trans[[ncol(trans) + 1]] <- NA
  
  trans <- data.frame(trans, transD, stringsAsFactors = FALSE)
  trans[[ncol(trans) + 1]] <- NA
  
  trans <- data.frame(trans, transE, stringsAsFactors = FALSE)
  trans <- cbind(trans, NA)
  
  names(trans) <- c("MutA Translation", "MutA Count", "MutA Proportion", "", "MutB (297 T->A)", "MutB Count", "MutB Proportion", "",
                       "MutC (368 T->N)", "MutC Count", "MutC Proportion", "", "MutD (382 G->R)", "Mut Count", "MutC Proportion", "",
                       "MutE", "MutE Count", "MutE Proportion", "")
  a4b7trans[[paste0(i)]]<-trans
  a4b7barcodes[[paste0(i)]] <- hold

}
names(a4b7trans) <- barcodesused
names(a4b7barcodes) <- barcodesused
a4b7Stats <- stats
write.csv(stats, "pathToA4B7Stats.csv", na = "", row.names = FALSE)

# Tuesday 2/5/19
# Adding in a version of V4 analysis that translates stuff to make Brandon happy

barcodesused <- a4b7primers

# Table to hold summary statistics
stats <- data.frame(Primer = character(), "# of 5' Indexing Sequences" = numeric(),
                    "Reads matching V4 reference" = numeric(), stringsAsFactors = FALSE)

#i = 1
V4translations <- list()
V4barcodes <- list()
for(i in a4b7primers){
  print(paste("Primer", i, barcode.list[i,2]))
  
  # Find the reads matching to this primer
  splitmatch=(barcodematches==i)&(!is.na(barcodematches))
  
  # pull the length of the read after the primer
  rev.seq.split.save <- str_sub(rev.seq[splitmatch], nchar(barcode.list[i,2]) + 1, nchar(rev.seq[splitmatch]))
  
  # Let's look for the reference sequence from Brandon, allowing 3 mismatches like we did in the app
  
  ref <- "CATATTAGACAAATAATCAACACTTGGCATAAAGTAGGCA"
  
  # Need to reverse complement the reverse strands
  rev.seq.fwd <- as.character(reverseComplement(DNAStringSet(rev.seq.split.save)))
  
  pos <- unlist(aregexec(ref, rev.seq.fwd, max.distance = 5))
  print(paste("Matches:", length(pos[pos > 0])))
  
  # We want the 114 bases upstream of this.
  v4 <- substr(rev.seq.fwd[pos > 0], pos[pos > 0] - 114, pos[pos > 0] - 1)
  
  # Make a table of the most common sequences
  table <- as.data.frame(table(v4))
  # Don't actually have this in codons yet, not sure what the correct frame is!
  table$AA <- as.character(translate(DNAStringSet(table$v4), if.fuzzy.codon = "solve"))
  if(length(table) == 2){
    V4translations[[i]] <- NA
    next
  }
  # And also a table of just the translations
  trans <- table[,2:3]
  trans <- aggregate(Freq ~ AA, trans, sum)
  
  # Discard singletons
  table <- table[(table$Freq > 1),]
  trans <- trans[(trans$Freq > 1),]
  
  # Calculate proportions
  table$proportion <- table$Freq/sum(table$Freq)
  trans$proportion <- trans$Freq/sum(trans$Freq)
  
  # Order counts by frequency
  table <- table[order(table$Freq, decreasing = TRUE),]
  trans <- trans[order(trans$Freq, decreasing = TRUE),]
  
  # Clean up tables a bit
  if(length(table[[1]]) > 0){
    table <- table[,c(1,3,2,4)]
    rownames(table) <- 1:length(table[[1]])
  }
  if(length(trans[[1]]) > 0){
    rownames(trans) <- 1:length(trans[[1]])
  }
  
  # Update our summary table
  stats[nrow(stats)+1,] <- c(paste0("A4B7.", barcode.list[i]), sum(splitmatch), length(pos[pos > 0]))
  
  # Need to make them the same length so we can output them in one file
  holdlen <- max(length(table[[1]]), length(trans[[1]]))
  
  # Actually making them one longer than the longest row so that nothing gets written over
  # (Fixing bug with Sirish method where the last row gets removed if the table is already the desired length)
  
  table[nrow(table) + (holdlen - nrow(table)) + 1,] <- NA
  trans[nrow(trans) + (holdlen - nrow(trans)) + 1,] <- NA
  
  # And now let's combine these into one dataframe
  table[[(ncol(table)) + 1]] <- " "
  hold <- data.frame(table, trans, stringsAsFactors = FALSE)
  
  colnames(hold) <- c("V4 Sequence", "V4 Translation", "V4 Count", "V4 Proportion", "", "V4 Translation", "V4 Count", "V4 Proportion")
  transHold <- hold[,6:8]
  V4translations[[as.character(i)]] <- transHold
  V4barcodes[[as.character(i)]] <- hold[,1:4]

  
}
names(V4translations) <- barcodesused
names(V4barcodes) <- barcodesused
colnames(stats) <- c("Primer", "# of 5' Indexing Sequences", "# of Reads Matching V4 Reference")
V4Stats <- stats
write.csv(stats, "pathToV4Stats.csv", na = "", row.names = FALSE)




## REFORMATTING BEGINS HERE -- If run info forms not properly set up - this is where the trickiness begins ##

## now want to create a wide dataframe per sample that shows reads/proportions 
force_bind = function(df1, df2) {
  colnames(df2) = colnames(df1)
  bind_rows(df1, df2)
}


## merge multiple dataframes into 1. One df for AA seq, one for sequences
fullDFs <- list()
fullnucDFs <- list()
# for each primer, run create a dataframe
for(i in V2V3primers){
  V2<- V2translations[[i]]
  V4 <- V4translations[[i]]
  a4b7 <- a4b7trans[[i]]
  if(ncol(a4b7 == 21)){
    a4b7 = a4b7[,1:20]
  }
  maxRows <- max(nrow(V2), nrow(V4), nrow(a4b7))
  if (is.data.frame(V2) && nrow(V2) < maxRows ){
  V2 <- force_bind(V2, as.data.frame(matrix(NA, nrow = maxRows-nrow(V2), ncol = 3)))
  }else if(!is.data.frame(V2)){V2 <- as.data.frame(matrix(nrow=maxRows, ncol = 3))}
  if (is.data.frame(V4) && nrow(V4) < maxRows){
  V4 <-force_bind(V4, as.data.frame(matrix(NA, nrow = maxRows-nrow(V4), ncol = 3)))
  }else if(!is.data.frame(V4)){V4 <- as.data.frame(matrix(nrow=maxRows, ncol = 3))}
  if(is.data.frame(a4b7) && nrow(a4b7) <maxRows){
  a4b7 <- force_bind(a4b7, as.data.frame(matrix(NA,nrow = maxRows-nrow(a4b7), ncol = 20)))
  }else if(!is.data.frame(a4b7)){a4b7 <- as.data.frame(matrix(nrow=maxRows, ncol = 3))}
  df <- cbind(V2, NA,V4,NA, a4b7)
  names(df) <- c("V2 Sequence", "V2 Count", "V2 Proportion", "", "V4 Translation", 
                 "V4 Count", "V4 Proportion", "", "MutA Translation", "MutA Count", 
                 "MutA Proportion", "", "MutB (295 T)", "MutB Count", "MutB Proportion", 
                 "", "MutC (368 T->N)", "MutC Count", "MutC Proportion", 
                 "", "MutD (382 G->R)", "MutD Count", "MutD Proportion", 
                 "", "Mut E", "MutE Count", "MutE Proportion")
  fullDFs[[i]] <- df
  
  V2nucs <- V2barcodes[[i]]
  V4nucs <- V4barcodes[[i]]
  a4b7nucs <- a4b7barcodes[[i]]
  
  maxRows <- max(nrow(V2nucs), nrow(V4nucs), nrow(a4b7nucs))
  if (is.data.frame(V2nucs) && nrow(V2nucs) < maxRows ){
    V2nucs <- force_bind(V2nucs, as.data.frame(matrix(NA, nrow = maxRows-nrow(V2nucs), ncol = 5)))
  }else if(!is.data.frame(V2nucs)){V2nucs <- as.data.frame(matrix(nrow=maxRows, ncol = 3))}
  if (is.data.frame(V4nucs) && nrow(V4nucs) < maxRows){
    V4nucs <-force_bind(V4nucs, as.data.frame(matrix(NA, nrow = maxRows-nrow(V4nucs), ncol = 4)))
  }else if(!is.data.frame(V4nucs)){V4nucs <- as.data.frame(matrix(nrow=maxRows, ncol = 3))}
  if(is.data.frame(a4b7nucs) && nrow(a4b7nucs) <maxRows){
    a4b7nucs <- force_bind(a4b7nucs, as.data.frame(matrix(NA,nrow = maxRows-nrow(a4b7nucs), ncol = 25)))
  }else if(!is.data.frame(a4b7nucs)){a4b7nucs <- as.data.frame(matrix(nrow=maxRows, ncol = 3))}
  fullnucDFs[[i]] <- cbind(V2nucs,NA, V4nucs, NA, a4b7nucs)
  
}
names(fullDFs) <- barcodesused



V2 <- list()
V4 <- list()
MutA <- list()
MutB <- list()
MutC <- list()
MutD <- list()
MutE <- list()
lapply(names(fullDFs), function(x){
  df <- fullDFs[[x]]
  V2[[x]] <<- df[,c(1:3)]
  V4[[x]] <<- df[,c(5:7)]
  MutA[[x]] <<- df[,c(9:11)]
  MutB[[x]] <<- df[,c(13:15)]
  MutC[[x]] <<- df[,c(17:19)]
  MutD[[x]] <<- df[,c(21:23)]
  MutE[[x]] <<- df[,c(25:27)]
})

concatDFs <- list()

for(i in names(fullDFs)){
  V2temp <- V2[[i]]
  V2temp <- V2temp[1:10,]
  V2temp <- V2temp[!is.na(V2temp$`V2 Sequence`),]
  V2temp$`V2 Proportion` <- V2temp$`V2 Count`/sum(V2temp$`V2 Count`)
  V4temp <- V4[[i]]
  V4temp <- V4temp[1:10,]
  V4temp <- V4temp[!is.na(V4temp$`V4 Translation`),]
  V4temp$`V4 Proportion` <- V4temp$`V4 Count`/sum(V4temp$`V4 Count`)
  MutAtemp <- MutA[[i]]
  MutAtemp$`MutA Count` <- as.numeric(MutAtemp$`MutA Count`)
  if(sum(is.na(MutAtemp$`MutA Translation`)) != nrow(MutAtemp)){
  MutAtemp <- aggregate(`MutA Count`~`MutA Translation`, MutAtemp, sum)
  MutAtemp <- MutAtemp[order(MutAtemp$`MutA Count`, decreasing = T),]
  MutAtemp$Proportion <- as.numeric(MutAtemp$`MutA Count`)/sum(as.numeric(MutAtemp$`MutA Count`))
  }
  if(nrow(MutAtemp) > 10){
    MutAtemp <- MutAtemp[1:10,]
    MutAtemp$Proportion <- as.numeric(MutAtemp$`MutA Count`)/sum(as.numeric(MutAtemp$`MutA Count`))
  }
  names(MutAtemp)[3] <- c("MutA Proportion")
  MutAtemp[,4] <- NA
  MutBtemp <- MutB[[i]]
  MutBtemp$`MutB Count` <- as.numeric(MutBtemp$`MutB Count`)
  if(sum(is.na(MutBtemp$`MutB (295 T)`)) != nrow(MutBtemp)){
  MutBtemp$`MutB (295 T)` <- substring(MutBtemp$`MutB (295 T)`, 4,6)
  MutBtemp <- aggregate(`MutB Count`~`MutB (295 T)`, MutBtemp, sum)
  MutBtemp <- MutBtemp[order(MutBtemp$`MutB Count`, decreasing = T),]
  MutBtemp$Proportion <- as.numeric(MutBtemp$`MutB Count`)/sum(as.numeric(MutBtemp$`MutB Count`))
  }
  if(nrow(MutBtemp) > 10){
    MutBtemp <- MutBtemp[1:10,]
    MutBtemp$Proportion <- as.numeric(MutBtemp$`MutB Count`)/sum(as.numeric(MutBtemp$`MutB Count`))
  }
  names(MutBtemp)[3] <- c("MutB Proportion")
  MutBtemp[,4] <- NA
  MutCtemp <- MutC[[i]]
  MutCtemp$`MutC Count` <- as.numeric(MutCtemp$`MutC Count`)
  if(sum(is.na(MutCtemp$`MutC (368 T->N)`)) != nrow(MutCtemp)){
  MutCtemp$`MutC (368 T->N)` <- substring(MutCtemp$`MutC (368 T->N)`, 6,8)
  MutCtemp <- aggregate(`MutC Count`~`MutC (368 T->N)`, MutCtemp, sum)
  MutCtemp <- MutCtemp[order(MutCtemp$`MutC Count`, decreasing = T),]
  MutCtemp$Proportion <- as.numeric(MutCtemp$`MutC Count`)/sum(as.numeric(MutCtemp$`MutC Count`))
  }
  if(nrow(MutCtemp) > 10){
    MutCtemp <- MutCtemp[1:10,]
    MutCtemp$Proportion <- as.numeric(MutCtemp$`MutC Count`)/sum(as.numeric(MutCtemp$`MutC Count`))
  }
  MutCtemp[,4] <- NA
  names(MutCtemp)[3] <- c("MutC Proportion")
  MutDtemp <- MutD[[i]]
  MutDtemp$`MutD Count` <- as.numeric(MutDtemp$`MutD Count`)
  MutDtemp$`MutD (382 G->R)` <- substring(MutDtemp$`MutD (382 G->R)`, 6,6)
  names(MutDtemp)[2] <- c("MutD Count")
  if(!is.na(MutDtemp$`MutD Count`[1])){
    MutDtemp <- aggregate(`MutD Count`~`MutD (382 G->R)`, MutDtemp, sum)
  }
  MutDtemp <- MutDtemp[order(MutDtemp$`MutD Count`, decreasing = T),]
  MutDtemp$`MutD Proportion` <- as.numeric(MutDtemp$`MutD Count`)/sum(as.numeric(MutDtemp$`MutD Count`))
  MutDtemp[,4] <- NA
  if(nrow(MutDtemp) > 10){
    MutDtemp <- MutDtemp[1:10,]
    MutDtemp$`MutD Proportion` <- as.numeric(MutDtemp$`MutD Count`)/sum(as.numeric(MutDtemp$`MutD Count`))
  }
  
  MutE1temp <- MutE[[i]]
  MutE1temp$`Mut E` <- substring(MutE1temp$`Mut E`, 6,6)
  if(!is.na(MutE1temp$`MutE Count`[1])){
  MutE1temp <- aggregate(`MutE Count`~`Mut E`, MutE1temp, sum)
  }
  MutE1temp$`MutE Count` <- as.numeric(MutE1temp$`MutE Count`)
  MutE1temp <- MutE1temp[order(MutE1temp$`MutE Count`, decreasing = T),]
  MutE1temp$`MutE Proportion` <- as.numeric(MutE1temp$`MutE Count`)/sum(as.numeric(MutE1temp$`MutE Count`))
  MutE1temp[,4] <- NA
  names(MutE1temp)[1:3] <- c("MutE1 Translation", "MutE1 Count", "MutE1 Proportion")
  if(nrow(MutE1temp) > 10){
    MutE1temp <- MutE1temp[1:10,]
    MutE1temp$`MutE1 Proportion` <- as.numeric(MutE1temp$`MutE1 Count`)/sum(as.numeric(MutE1temp$`MutE1 Count`))
  }
  MutE2temp <- MutE[[i]]
  MutE2temp$`Mut E` <- substring(MutE2temp$`Mut E`, 8,8)
  if(!is.na(MutE2temp$`MutE Count`[1])){
    MutE2temp <- aggregate(`MutE Count`~`Mut E`, MutE2temp, sum)
  }
  MutE2temp$`MutE Count` <- as.numeric(MutE2temp$`MutE Count`)
  MutE2temp <- MutE2temp[order(MutE2temp$`MutE Count`, decreasing = T),]
  MutE2temp$`MutE Proportion` <- as.numeric(MutE2temp$`MutE Count`)/sum(as.numeric(MutE2temp$`MutE Count`))
  MutE2temp[,4] <- NA
  names(MutE2temp)[1:3] <- c("MutE2 Translation", "MutE2 Count", "MutE2 Proportion")
  if(nrow(MutE2temp) > 10){
    MutE2temp <- MutE2temp[1:10,]
    MutE2temp$`MutE2 Proportion` <- as.numeric(MutE2temp$`MutE2 Count`)/sum(as.numeric(MutE2temp$`MutE2 Count`))
  }
  concatDFs[[i]] <- qpcR:::cbind.na(V2temp, V4temp, MutAtemp, MutBtemp, MutCtemp, MutDtemp, MutE1temp, MutE2temp)
  names(concatDFs[[i]]) <- c("V2 Sequence", "V2 Count", "V2 Proportion", "V4 Translation", 
                             "V4 Count", "V4 Proportion", "MutA Translation", "MutA Count", 
                             "MutA Proportion", "", "MutB (295 N-linked glycan)", "MutB Count", "MutB Proportion", 
                             "", "MutC (371 N-linked glycan)", "MutC Count", "MutC Proportion", 
                             "", "MutD (382 G->R)", "MutD Count", "MutD Proportion", 
                             "", "Mut E1", "MutE1 Count", "MutE1 Proportion",
                             "", "Mut E2", "MutE2 Count", "MutE2 Proportion", "")
}
# construct empty DFs with proper column names
V2 <- concatDFs[[1]][0,1:3]
V4 <- concatDFs[[1]][0,4:6]
MutA <- concatDFs[[1]][0,7:9]
MutB <- concatDFs[[1]][0,11:13]
MutC <- concatDFs[[1]][0,15:17]
MutD <- concatDFs[[1]][0,19:21]
MutE1 <- concatDFs[[1]][0,23:25]
MutE2 <- concatDFs[[1]][0,27:29]

## if this fails, then you have a different number of primers than you do names for some reason
names(concatDFs) <- names

lapply(names, function(x){
  V2 <<- rbind(V2, cbind(concatDFs[[x]][,1:3], Animal=rep(x, nrow(concatDFs[[x]]))))
  V4 <<- rbind(V4, cbind(concatDFs[[x]][,4:6], Animal=rep(x, nrow(concatDFs[[x]]))))
  MutA <<- rbind(MutA, cbind(concatDFs[[x]][,7:9], Animal=rep(x, nrow(concatDFs[[x]]))))
  MutB <<- rbind(MutB, cbind(concatDFs[[x]][,11:13], Animal=rep(x, nrow(concatDFs[[x]]))))
  MutC <<- rbind(MutC, cbind(concatDFs[[x]][,15:17], Animal=rep(x, nrow(concatDFs[[x]]))))
  MutD <<- rbind(MutD, cbind(concatDFs[[x]][,19:21], Animal=rep(x, nrow(concatDFs[[x]]))))
  MutE1 <<- rbind(MutE1, cbind(concatDFs[[x]][,23:25], Animal=rep(x, nrow(concatDFs[[x]]))))
  MutE2 <<- rbind(MutE2, cbind(concatDFs[[x]][,27:29], Animal=rep(x, nrow(concatDFs[[x]]))))
})

DFlist <- list(V2,V4,MutA,MutB,MutC,MutD,MutE1,MutE2)
DFnames <- c("V2", "V4", "MutA", "MutB", "MutC","MutD", "MutE1", "MutE2")
names(DFlist) <- DFnames
listMatrix <- lapply(DFlist, function(x){
  x <<- x
  x <- x[complete.cases(x),]
  if(nrow(x) == 0){
    return(NA)
  }
  uniqueAA <- unique(x[,1])
  uniqueAnimal <- unique(x[,4])
  AAdf <- as.data.frame(uniqueAA)
  for (curAnimal in names){
    AAdf[,paste(curAnimal)] <- NA
    for (i in 1:nrow(AAdf)){
      curAA <- AAdf$uniqueAA[i]
      matchIndex <- which(x[,1] == curAA & x[,4] == curAnimal)
      if (length(matchIndex) > 0){
        AAdf[i,paste(curAnimal)] <- x[,3][matchIndex] 
      }
    }
  }
  AAdf[is.na(AAdf)] <- 0
  sums <- rowSums(AAdf[,2:(length(names)+1)], na.rm=TRUE)
  AAdf <- AAdf[order(sums, decreasing =T),]
  return(AAdf)
})
names(listMatrix) <- DFnames
if(sum(is.na(listMatrix))> 0){
  listMatrix <- listMatrix[-is.na(listMatrix)]
}

runInfoList <- merge(runInfoList, V2Stats, by.x = "Barcodes", by.y = "Primer", all.x = T)
runInfoList <- merge(runInfoList, a4b7Stats, by.x = "Barcodes", by.y = "Primer", all.x = T)
runInfoList <- merge(runInfoList, V4Stats, by.x = "Barcodes", by.y = "Primer", all.x = T)
wb<- createWorkbook()
addWorksheet(wb, "Run Info Form")
writeData(wb, "Run Info Form", runInfoList)
lapply(names(listMatrix), function(x){
  addWorksheet(wb, x)
  writeData(wb, x, listMatrix[[x]])
})


lapply(names(concatDFs), function(dfName){
  if(nchar(dfName) > 31){
    wsName <- substring(dfName, 1,31)
  }else{
    wsName <- dfName
  }
  addWorksheet(wb, wsName)
  writeData(wb, sheet = wsName, concatDFs[[dfName]])
  setColWidths(wb, wsName, cols = 1:ncol(concatDFs[[dfName]]), widths = "auto")
})

saveWorkbook(wb, "pathToFinalExcelOutputV2V3A4B7.xlsx", overwrite = T)
