require(plyr)

KOC_SAMPLES <- c('N406', 'N411', 'N416') #SWP
OC_SAMPLES <- c('N407', 'N412', 'N417') #OCWP
KCC_SAMPLES <- c('N406', 'N411', 'N416', 'N409', 'N414', 'N419') #SWP SMP
KCCW_SAMPLES <- c('N406', 'N411', 'N416') #SWP
KCCM_SAMPLES <- c('N409', 'N414', 'N419') #SMP
CC_SAMPLES <- c('N408', 'N410', 'N415', 'N420' , 'N418','N413') #CCWP CCMP
CCW_SAMPLES <- c('N408', 'N418','N413') #CCWP
CCM_SAMPLES <- c('N410', 'N415', 'N420') #CCMP
ALL_SAMPLES <- unique(c(KOC_SAMPLES, KCC_SAMPLES, OC_SAMPLES, CC_SAMPLES))
NORM_SWATH <- FALSE
LOG_BASE <- 2
FC_TH <- 1

PVREQ <- 0.05

tst <- function(row, k, e)
{
    a = length(na.omit(row[c(k)]))
    b = length(na.omit(row[c(e)]))
    t = list(p.value = NA) # If test cannot be performed
    try(
{
    t <- wilcox.test(as.numeric(na.omit(row[c(k)])),
                     as.numeric(na.omit(row[c(e)])),
                     correct = F, exact = F)
}, TRUE)
if(is.finite(t$p.value))
{
    return(t$p.value)
}
else
{
    return(1)
}
}
##############
swath_rt_data <- read.table("data/ions.txt", header = T, sep = '\t', stringsAsFactors = FALSE)
swath_frg_data <- data.frame(Pep = swath_rt_data$Peptide, 
                            RT = swath_rt_data$RT,
                            precursor_charge = swath_rt_data$Precursor.Charge, 
                            fragment_series = swath_rt_data$Ion.Type,
                            fragment_number = swath_rt_data$Residue,
                            fragment_charge = swath_rt_data$Fragment.Charge,
                            PrecursorMZ = swath_rt_data$Precursor.MZ,
                            FragmentMZ = swath_rt_data$Fragment.MZ,
                            stringsAsFactors = FALSE)
swath_frg_data$ID = paste(sep = "_", swath_frg_data$Pep, swath_frg_data$precursor_charge, 
                           swath_frg_data$fragment_series, swath_frg_data$fragment_number, 
                           swath_frg_data$fragment_charge)
frgs <- read.csv("fragments.txt", sep="\t", header = T, stringsAsFactors = F)
frgs <- rename(frgs, c("clean_peptide_sequence" = "Pep" ))
frgs$ID = paste(sep="_", frgs$Pep, frgs$precursor_charge, 
                           frgs$fragment_series, frgs$fragment_number, 
                           frgs$fragment_charge)
frgs <- unique(frgs)
swath_frg_data <- unique(swath_frg_data)
frgs_ <- merge(swath_frg_data, frgs, by = "ID", all.y=T, all.x=F)
frgs_ <- na.omit(frgs_)
frgs_ <- unique(frgs_)
frgs_$Pep <- frgs_$Pep.y
frgs_ <- frgs_[,c("Pep", "RT", "PrecursorMZ", "FragmentMZ")]
write.table(frgs_,"fragment_base.txt", quote = F, sep = "\t", row.names=F)
swath_rt_data <- data.frame(Pep = swath_rt_data$Peptide, 
                            RT = swath_rt_data$RT, Prot = swath_rt_data$Protein, 
                            stringsAsFactors = FALSE)
swath_rt_data$Pep <- gsub(x= swath_rt_data$Pep, pattern= '\\[.*?\\]', replacement= '')
swath_rt_data$Pep <- gsub(x= swath_rt_data$Pep, pattern= '-', replacement= '', fixed= TRUE)

swath_prot <- swath_rt_data[,c("Pep", "Prot")]
swath_prot <- swath_prot[!duplicated(swath_prot$Pep),]
swath_prot$id <- sapply(strsplit(swath_prot$Prot, "|", fixed = TRUE), "[", 2)

#tmp_p <- swath_prot[,c("Prot", "id")]
#tmp_p <- unique(tmp_p)

maps <- read.csv("./genetab.tab", sep="\t", header = T)
maps <- maps[, c("Entry", "Protein.names")]
maps <- rename(maps, c("Entry" = "id", "Protein.names" = "Description"))
swath_prot <- merge(swath_prot, maps, all.x = TRUE)


#require(UniProt.ws)
#taxId(UniProt.ws) <- 9606

#sl <- select(UniProt.ws, tmp_p$id[1:3], "ENSEMBL", "UNIPROTKB")

#require(biomaRt)
#uniprotV <- useMart("unimart", dataset="uniprot")
#d <- getBM(attributes=c("accession", "name"), mart = uniprotV)


swath_rt_data <- aggregate(RT ~ Pep, data = swath_rt_data, FUN = mean)
swath_rt_data <- merge(swath_rt_data, swath_prot, by = "Pep", all=T)
swath_rt_data <- unique(swath_rt_data)

##############

peps <- read.table("data//peptides1_wo407.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
#peps <- read.table("data//peptides_0_rec.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
#pepsx = merge(peps, peps0, by="clean_peptide_sequence")

if (NORM_SWATH)
{
    p <- peps[,-1]
    avgs <- apply(p, 2, mean)
    p <- sweep(p, 2, avgs, "/") 
    p$Pep <- peps$clean_peptide_sequence
    peps <- p
    rm(p)
} else
{
    peps <- rename(peps, c('clean_peptide_sequence' = 'Pep'))    
}

scaf_data <- read.table("Peptide Report for 2_loaded.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
scaf_ok <- data.frame(Pep = scaf_data$Peptide.sequence, stringsAsFactors = F)
scaf_ok <- unique(scaf_ok)
scaf_okl <- scaf_ok$Pep

peps$Scaf <- sapply(peps$Pep, FUN = function(x) x %in% scaf_okl)

#require(biomaRt)
#mart = useMart("ensembl")
#mart <- useMart("unimart", dataset="uniprot")
#bm <- getBM(attributes=c("accession", "name"), mart = mart)
#library(UniProt.ws)
#taxId(UniProt.ws) <- 9606



####
peps$CCW <- apply(peps[,CCW_SAMPLES], 1, FUN = mean)
peps$KCCW <- apply(peps[,KCCW_SAMPLES], 1, FUN = mean)
peps$pepsCCWtoKCCW <- peps$CCW/peps$KCCW
peps$ccwlfc <- log(peps$pepsCCWtoKCCW, base = LOG_BASE)
peps$ccwpv <- apply(peps, 1, tst, KCCW_SAMPLES, CCW_SAMPLES) 
peps$ccwqv <- p.adjust(peps$ccwpv, method = "BH")
ccwlist <- unique( subset(peps, ccwpv < PVREQ & abs(ccwlfc) > 1, c(Pep, ccwlfc, ccwpv, Scaf)))
rownames(ccwlist) <- NULL
ccwlist <- merge(ccwlist, swath_rt_data, by = "Pep", all.x =TRUE)
ccwlist <- ccwlist[order(-abs(ccwlist$ccwlfc), ccwlist$ccwpv),]
ccwlist <- rename(ccwlist, c("ccwlfc" = "log2fc", "ccwpv" = "pv", "Prot" = "Protein"))
write.table(ccwlist[,c("Pep", "RT", "log2fc", "pv", "Protein", "Description")], 
            "CCW_peptides.txt", sep = '\t', row.names = F, quote = F)
####
peps$CCM <- apply(peps[,CCM_SAMPLES], 1, FUN = mean)
peps$KCCM <- apply(peps[,KCCM_SAMPLES], 1, FUN = mean)
peps$pepsCCMtoKCCM <- peps$CCM/peps$KCCM
peps$ccmlfc <- log(peps$pepsCCMtoKCCM, base = LOG_BASE)
peps$ccmpv <- apply(peps, 1, tst, KCCM_SAMPLES, CCM_SAMPLES) 
peps$ccmqv <- p.adjust(peps$ccmpv, method = "BH")
ccmlist <- unique( subset(peps, ccmpv < PVREQ & abs(ccmlfc) > 1, c(Pep, ccmlfc, ccmpv, Scaf)))
rownames(ccmlist) <- NULL
ccmlist <- merge(ccmlist, swath_rt_data, by = "Pep", all.x =TRUE)
ccmlist <- ccmlist[order(-abs(ccmlist$ccmlfc), ccmlist$ccmpv),]
ccmlist <- rename(ccmlist, c("ccmlfc" = "log2fc", "ccmpv" = "pv", "Prot" = "Protein"))
write.table(ccmlist[,c("Pep", "RT", "log2fc", "pv", "Protein", "Description")], 
            "CCM_peptides.txt", sep = '\t', row.names = F, quote = F)


####
peps$CC <- apply(peps[,CC_SAMPLES], 1, FUN = mean)
peps$KCC <- apply(peps[,KCC_SAMPLES], 1, FUN = mean)
peps$pepsCCtoKCC <- peps$CC/peps$KCC
peps$cclfc <- log(peps$pepsCCtoKCC, base = LOG_BASE)
peps$ccpv <- apply(peps, 1, tst, KCC_SAMPLES, CC_SAMPLES) 
peps$ccqv <- p.adjust(peps$ccpv, method = "BH")
cclist <- unique( subset(peps, ccpv < PVREQ & abs(cclfc) > 1, c(Pep, cclfc, ccpv, Scaf)))
rownames(cclist) <- NULL
cclist <- merge(cclist, swath_rt_data, by = "Pep", all.x =TRUE)
cclist <- cclist[order(-abs(cclist$cclfc), cclist$ccpv),]
cclist <- rename(cclist, c("cclfc" = "log2fc", "ccpv" = "pv", "Prot" = "Protein"))
cclist$M <- sapply(cclist$Pep, FUN = function(x) x %in% ccmlist$Pep)
cclist$W <- sapply(cclist$Pep, FUN = function(x) x %in% ccwlist$Pep)
write.table(cclist[,c("Pep", "RT", "log2fc", "pv", "Protein", "W", "M", "Description")], 
            "CC_peptides.txt", sep = '\t', row.names = F, quote = F)

ccmlist <- ccmlist[!(ccmlist$Pep %in% cclist$Pep),]
write.table(ccmlist[,c("Pep", "RT", "log2fc", "pv", "Protein", "Description")], 
            "CCMuniq_peptides.txt", sep = '\t', row.names = F, quote = F)
ccmlist$M <- TRUE
ccmlist$W <- FALSE
ccmtlist <- rbind(cclist, ccmlist)
ccmtlist <- ccmtlist
write.table(ccmtlist[,c("Pep", "RT", "log2fc", "pv", "Protein", "Description")], 
            "CCMtotal.txt", sep = '\t', row.names = F, quote = F)

ccwlist <- ccwlist[!(ccwlist$Pep %in% cclist$Pep),]
write.table(ccwlist[,c("Pep", "RT", "log2fc", "pv", "Protein", "Description")], 
            "CCWuniq_peptides.txt", sep = '\t', row.names = F, quote = F)
ccwlist$M <- FALSE
ccwlist$W <- TRUE
ccwtlist <- rbind(cclist, ccwlist)
ccwtlist <- ccmtlist
write.table(ccwtlist[,c("Pep", "RT", "log2fc", "pv", "Protein", "Description")], 
            "CCWtotal.txt", sep = '\t', row.names = F, quote = F)
####
peps$OC <- apply(peps[,OC_SAMPLES], 1, FUN = mean)
peps$KOC <- apply(peps[,KOC_SAMPLES], 1, FUN = mean)
peps$pepsOCtoKOC <- peps$OC/peps$KOC
peps$oclfc <- log(peps$pepsOCtoKOC, base = LOG_BASE)
peps$ocpv <- apply(peps, 1, tst, KOC_SAMPLES, OC_SAMPLES) 
peps$ocqv <- p.adjust(peps$ocpv, method = "BH")
oclist <- unique(subset(peps, ocpv < PVREQ & abs(oclfc) > 1, c(Pep, oclfc, ocpv)))
rownames(oclist) <- NULL
oclist <- merge(oclist, swath_rt_data, by = "Pep", all.x =TRUE)
oclist <- oclist[order(-abs(oclist$oclfc), oclist$ocpv),]
oclist <- rename(oclist, c("oclfc" = "log2fc", "ocpv" = "pv", "Prot" = "Protein"))
write.table(oclist[,c("Pep", "RT", "log2fc", "pv", "Protein", "Description")], 
            "OC_peptides.txt", sep = '\t', row.names = F, quote = F)
###########################


both <- intersect(oclist$Pep, cclist$Pep)
difr <- setdiff(oclist$Pep, cclist$Pep)
AA_SAMPLES <- c(CC_SAMPLES, OC_SAMPLES)
KAA_SAMPLES <- c(KCC_SAMPLES, KOC_SAMPLES)

peps$AA <- apply(peps[,AA_SAMPLES], 1, FUN = mean)
peps$KAA <- apply(peps[,KAA_SAMPLES], 1, FUN = mean)
peps$pepsAAtoKAA <- peps$AA/peps$KAA
peps$AAlfc <- log(peps$pepsAAtoKAA, base = LOG_BASE)
peps$AApv <- apply(peps, 1, tst, KAA_SAMPLES, AA_SAMPLES) 
peps$AAqv <- p.adjust(peps$AApv, method = "BH")
AAlist <- unique(subset(peps, AApv < PVREQ & abs(AAlfc) > 1, c(Pep, AAlfc, AApv)))
rownames(AAlist) <- NULL
AAlist <- AAlist[order(-abs(AAlist$AAlfc), AAlist$AApv),]
AAlist <- merge(AAlist, swath_rt_data, by = "Pep", all.x =TRUE)
write.table(AAlist[,c("Pep", "RT")], "AA_peptides.txt", sep = '\t', row.names = F)


pconst <- peps[(abs(peps$cclfc) < 0.1  & abs(peps$oclfc) < 0.1),]
#pconst <- pconst[c(10,12),]


pconstg <- peps
rownames(pconstg) <- peps$Pep
pconstg <- pconstg[,ALL_SAMPLES]
avgs <- apply(pconstg, 1, mean)
pconstg <- sweep(pconstg, 1, avgs, "-")
pconstg <- sweep(pconstg, 1, avgs, "/")
pconstg <- pconstg^2
sqs <- apply(pconstg, 1, sum)
sqs <- sqs[order(sqs)]

constp <- names(sqs[1:30])
constp <- peps[peps$Pep %in% constp,]
constp <- constp[constp$ocpv > 0.05 & constp$ccpv > 0.05 & constp$ccwpv > 0.05 & constp$ccmpv > 0.05,]
constp <- constp[constp$Scaf,]
for (i in 1:nrow(constp))
{
    barplot(unlist(constp[i,ALL_SAMPLES]))
}
constp <- merge(constp, swath_rt_data, by = "Pep", all.x =TRUE)
#constp <- merge(constp, swath_prot, by = "Pep", all.x =TRUE)

constp <- rename(constp, c("Prot" = "Protein"))

write.table(constp[,c("Pep", "RT", "Protein")], 
            "const_peptides.txt", sep = '\t', row.names = F, quote = F)



#############################
stop("Run ended.")

peps1 <- read.table("data//peptides_rec_wo407.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
if (NORM_SWATH)
{
    p <- peps1[,-1]
    avgs <- apply(p, 2, mean)
    p <- sweep(p, 2, avgs, "/") 
    p$Pep <- peps1$clean_peptide_sequence
    peps1 <- p
    rm(p)
} else
{
    peps1 <- rename(peps1, c('clean_peptide_sequence' = 'Pep'))    
}


peps1$CC <- apply(peps1[,CC_SAMPLES], 1, FUN = mean)
peps1$KCC <- apply(peps1[,KCC_SAMPLES], 1, FUN = mean)
peps1$peps1CCtoKCC <- peps1$CC/peps1$KCC
peps1$cclfc <- log(peps1$peps1CCtoKCC, base = LOG_BASE)
peps1$ccpv <- apply(peps1, 1, tst, KCC_SAMPLES, CC_SAMPLES) 
peps1$ccqv <- p.adjust(peps1$ccpv, method = "BH")
cclist1 <- unique( subset(peps1, ccpv < PVREQ & abs(cclfc) > 1, c(Pep, cclfc, ccpv)))
rownames(cclist1) <- NULL
cclist1 <- cclist1[order(-abs(cclist1$cclfc), cclist1$ccpv),]
cclist1 <- merge(cclist1, swath_rt_data, by = "Pep", all.x =TRUE)
write.table(cclist1[,c("Pep", "RT")], "CC_peptides1.txt", sep = '\t', row.names = F)

peps1$OC <- apply(peps1[,OC_SAMPLES], 1, FUN = mean)
peps1$KOC <- apply(peps1[,KOC_SAMPLES], 1, FUN = mean)
peps1$peps1OCtoKOC <- peps1$OC/peps1$KOC
peps1$oclfc <- log(peps1$peps1OCtoKOC, base = LOG_BASE)
peps1$ocpv <- apply(peps1, 1, tst, KOC_SAMPLES, OC_SAMPLES) 
peps1$ocqv <- p.adjust(peps1$ocpv, method = "BH")
oclist1 <- unique(subset(peps1, ocpv < PVREQ & abs(oclfc) > 1, c(Pep, oclfc, ocpv)))
rownames(oclist1) <- NULL
oclist1 <- oclist1[order(-abs(oclist1$oclfc), oclist1$ocpv),]
oclist1 <- merge(oclist1, swath_rt_data, by = "Pep", all.x =TRUE)
write.table(oclist1[,c("Pep", "RT")], "OC_peptides1.txt", sep = '\t', row.names = F)


both1 <- intersect(oclist1$Pep, cclist1$Pep)
difr1 <- setdiff(oclist1$Pep, cclist1$Pep)