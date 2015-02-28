rm(list=ls(all=TRUE))

set.seed(3)
fake1<-matrix(c(rbinom(10, 1, 0.5),rbinom(10, 1, 0.8), rbinom(10, 1, 0.9),
         rbinom(10, 1, 0.1), rbinom(10, 1, 0.2), rbinom(10, 1, 0.8), 
         rbinom(10, 1, 0.4)), nrow=7, dimnames=list(c("a1","a2","a3",
                                                      "a4","a5","a6","a7"),
                                                    c("sp1","sp2","sp3","sp4",
                                                      "sp5","sp6","sp7","sp8",
                                                      "sp9","sp10")))

library(vegan)
library(magrittr)
library(dplyr)

fake1
rowSums(fake1)

#beta diversity based on pairwise differences in species richness
beta.rich <- function(.dados, .reps){
  reps<-.reps
  jac0<-as.data.frame(matrix(0, ncol=nrow(.dados), nrow=nrow(.dados)), row.names = row.names(.dados))
  jac1<-as.data.frame(matrix(0, ncol=nrow(.dados), nrow=nrow(.dados)), row.names = row.names(.dados))
  for(k in 1:reps) {
    for(j in 1:nrow(.dados)) {
      for(i in j:nrow(.dados)) {
        if (sum(.dados[i,])==sum(.dados[j,])) {
          r1<-as.numeric(vegdist(.dados[c(i,j),], method="jac"))
          jac0[i,j]<-r1
        } else if (sum(.dados[i,]) > sum(.dados[j,])) {
          trans1<-rrarefy(.dados[i,], sample = sum(.dados[j,]))
          r2<-as.numeric(vegdist(rbind(trans1, .dados[j,]), method="jac"))
          jac0[i,j]<-r2
        } else if (sum(.dados[i,]) < sum(.dados[j,])) {
          trans2<-rrarefy(.dados[j,], sample = sum(.dados[i,]))
          r3<-as.numeric(vegdist(rbind(.dados[i,],trans2), method="jac"))
          jac0[i,j]<-r3
        }
      }
    }
    jac1<-jac1+jac0
  }
  as.dist(jac1/.reps)
}

beta.rich(fake1, 999)

#nrow=(nrow(.dados)*(((nrow(.dados)^2)-nrow(.dados))/2)*(.reps))
#this is the number of elements in a vector made from a distance matrix
#multiplied by the number of times you repeat the procedure

# extract summary statistics from randomization ---------------------------

renamer <- function(x, pattern, replace) {
  for (i in seq_along(pattern))
    x <- gsub(pattern[i], replace[i], x)
  x
}

summary.beta <- function(.dados, .reps){
  reps<-.reps
  jac0<-as.data.frame(matrix(0, ncol=nrow(.dados), nrow=nrow(.dados)), row.names = row.names(.dados), colnames = row.names(.dados))
  jac2<-matrix(0, ncol=3, dimnames = list(c(), c("site1", "site2", "betadiver")))
  for(k in 1:reps) {
    jac1 <- matrix(0, nrow = nrow(.dados), ncol = nrow(.dados)) %>%
      lower.tri %>%
      which(arr.ind = TRUE)
    jac1<-as.data.frame(jac1)
    jac1$betadiver<-rep(0, nrow(jac1))
    jac1<-jac1 %>%
      mutate(site1 = renamer(row, seq(1:nrow(.dados)), row.names(.dados)), 
             site2 = renamer(col, seq(1:nrow(.dados)), row.names(.dados))) %>%
      select(-row, -col)
    for(j in 1:nrow(.dados)) {
      for(i in j:nrow(.dados)) {
        if (sum(.dados[i,])==sum(.dados[j,])) {
          r1<-as.numeric(vegdist(.dados[c(i,j),], method="jac"))
          jac0[i,j]<-r1
        } else if (sum(.dados[i,]) > sum(.dados[j,])) {
          trans1<-rrarefy(.dados[i,], sample = sum(.dados[j,]))
          r2<-as.numeric(vegdist(rbind(trans1, .dados[j,]), method="jac"))
          jac0[i,j]<-r2
        } else if (sum(.dados[i,]) < sum(.dados[j,])) {
          trans2<-rrarefy(.dados[j,], sample = sum(.dados[i,]))
          r3<-as.numeric(vegdist(rbind(.dados[i,],trans2), method="jac"))
          jac0[i,j]<-r3
        }
      }
    }
    jac1$betadiver <- as.vector(as.dist(jac0))
    jac2 <- rbind(jac2, jac1)         
  }
  jac2[-1,]
}


#rarefy communities to the same levels of species richness and then run jaccard dissimilarity
beta.rare<-function(dados, reps, comp) {
  reps<-reps
  jac0<-as.data.frame(matrix(0, ncol=nrow(dados), nrow=nrow(dados)))
  for(i in 1:reps) {
  r1<-rrarefy(dados, sample = comp)
  jac1<-as.data.frame(as.matrix(vegdist(r1, method="jac")))
  jac0<-jac0+jac1
}
  r2<-jac0/reps
  rownames(r2)<-rownames(dados)
  colnames(r2)<-rownames(dados)
  as.dist(r2)
}