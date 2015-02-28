#this function calculates the dissimilarity between pairs of sites, but taking into account differences
#in species richness in each comparison
#in order to improve the estimate of the rarefied dissimilarity, this rarefaction should be run for a 
#number of times, from which we take the average dissimilarity between pairs of communities after
#.reps runs
beta.rich <- function(.dados, .reps){
  reps<-.reps
  #create an empty data frame where you are going to store each pairwise dissimilarity value for each run
  jac0<-as.data.frame(matrix(0, ncol=nrow(.dados), nrow=nrow(.dados)), row.names = row.names(.dados)) 
  #create an empty data frame where you are going to sum up the values for each pairwise comparison
  jac1<-as.data.frame(matrix(0, ncol=nrow(.dados), nrow=nrow(.dados)), row.names = row.names(.dados))
  for(k in 1:reps) {
    for(j in 1:nrow(.dados)) {
      for(i in j:nrow(.dados)) {
        #for sites that match their species richness, there is no need of rarefation in the comparison
        if (sum(.dados[i,])==sum(.dados[j,])) { 
          r1<-as.numeric(vegdist(.dados[c(i,j),], method="jac"))
          jac0[i,j]<-r1
        } else if (sum(.dados[i,]) > sum(.dados[j,])) { #if one sites ir richer than the other, apply rarefaction
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
    jac1<-jac1+jac0 #for each run, sum the values obtained in one matrix
  }
  as.dist(jac1/.reps) #take the average for the number of runs
}