library(ape)
ins1<- read.GenBank(c("DQ250563.1", "J00747.1", "NM_008386.4",
               "XM_032892182.1", "DQ250565.1", "XM_031385017.1",
               "XM_021204825.1", "XM_036176265.1",
               "XM_013112606.3", "AJ009655.1"))

write.dna( ins1, file = 'ins1.fasta', format = 'fasta' )

coi <- read.GenBank(c("JF445255.1", "MZ661173.1", "JF499320.1", 
                      "MZ364426.1", "KY605336.1", "JQ667769.1", 
                      "JF445031.1", "HQ980049.1", "JF444326.1",
                      "KC750830.1"))

write.dna(coi, file = 'coi.fasta', format = 'fasta' )

library(DECIPHER)
all_ins1 <- "ins1.fasta"
all_coi <- "coi.fasta"

dna_all_ins1 <- readDNAStringSet(all_ins1)
dna_all_coi <- readDNAStringSet(all_coi)

ins1_align <- AlignSeqs(dna_all_ins1)
coi_align <- AlignSeqs(dna_all_coi)


BrowseSeqs(ins1_align)
BrowseSeqs(coi_align)


writeXStringSet(ins1_align, file="ins1_aligned.fasta")
writeXStringSet(coi_align, file="coi_aligned.fasta")


#NJ



library(phangorn)

all.phy <- read.phyDat('all_alignment.fasta', format = 'fasta', type = 'DNA')

all.phy
all.align <- read.dna('all_alignment.fasta', format = 'fasta')

dist.matrix <- dist.dna(all.align)
dist.matrix


tree.guide <- NJ(dist.matrix)
mod.test <- modelTest(all.phy, tree = tree.guide)
mod.test

# I used AIC to choose a model. My choice was SYM+G(4)+I but i could use it further 
# that is why i chose F81

dist.matrix <- dist.dna(all.align, model = "F81")
tree <- NJ(dist.matrix)
tree
tree.root <- root(tree, outgroup = 'Homo sapiens')

plot(tree.root, type = "phylogram", main = 'Neighbor Joining, INS-1 & COI')

write.tree(tree.root, file = 'nj_tree.tre')




# Parsimony
library(phylotools)

bab.tree <- bab(all.align)
bab.root <- root(bab.tree, 
                 outgroup = c('Homo sapiens'))
plot(bab.tree, main = 'Parsimony, BAB')

tree.SPR <- pratchet(all.align, maxit = 10000, minit = 100, k = 10,
                     all = T, rearrangements = 'SPR', trace = 0)

tree.NNI <- pratchet(all.align, maxit = 10000, minit = 100, k = 10,
                     all = T, rearrangements = 'NNI', trace = 0)
tree.SPR <- acctran(tree.SPR, all.phy)
tree.NNI <- acctran(tree.NNI, all.phy)

SPR.root <- root(tree.SPR, 
                 outgroup = c('Homo sapiens'))
NNI.root <- root(tree.NNI, 
                 outgroup = c('Homo sapiens'))

plot(SPR.root, main = 'Parsimony, SPR')
plot(NNI.root, main = 'Parsimony, NNI')

parsimony(bab.root, all.phy)
parsimony(SPR.root, all.phy)
parsimony(NNI.root, all.phy)

plotBS(SPR.root[[2]], type = "p", main = 'Parsimony, SPR')
plotBS(NNI.root[[2]], type = "p", main = 'Parsimony, SPR')




# Maximum Likelihood
dist <- dist.ml(all.align)
nj.tree <- nj(dist)

fit <- pml(nj.tree,data = all.phy)
fit

fitJC <- optim.pml(fit, rearrangement = "NNI")
plot(fitJC$tree, main = "JC, NNI rearrangement")



fitGTR.G <- update(fit, model = "GTR", k = 4)
fitGTR.G <- optim.pml(fitGTR.G, model = "GTR", optGamma = T, rearrangement = "stochastic", control = pml.control(trace = 0))
bs <- bootstrap.pml(fitGTR.G, bs=100, optNni=TRUE, control = pml.control(trace = 0))

tree.root <- root(fitGTR.G$tree, outgroup = c('Homo sapiens'))

plotBS(tree.root, bs, main = "GTR model. Bootstrap", type = "p",
       bs.col="red", p = 0.5, digits = 2)
