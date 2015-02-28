#run this function so that it renames row and columns of your lower.tri to match that
#of your original community matrix
renamer <- function(x, pattern, replace) {
  for (i in seq_along(pattern))
    x <- gsub(pattern[i], replace[i], x)
  x
}

#this function calculates the dissimilarity between pairs of sites, but taking into account differences
#in species richness in each comparison, n times
#this function stores the result of each randomization in a data frame, so that you can do some summary
#statistics with that distribution per comparison, or summarise it all and throw it back into a distance
#matrix
summary.beta <- function(.dados, .reps){
  reps<-.reps
  #create a dataframe to store the information in each run
  jac0<-as.data.frame(matrix(0, ncol=nrow(.dados), nrow=nrow(.dados)))
  #create a base dataframe where you are going to bind all the results at the end
  jac2<-matrix(0, ncol=3, dimnames = list(c(), c("site1", "site2", "betadiver")))
  for(k in 1:reps) {
    #create a data frame where you are going to store the pairwise results for each run
    jac1 <- matrix(0, nrow = nrow(.dados), ncol = nrow(.dados)) %>%
      lower.tri %>% #get only the lower triangle of the distance matrix
      which(arr.ind = TRUE) #save that triangle as a matrix with two columns representing the pairwise comparison
    jac1<-as.data.frame(jac1) #turn that into a data frame
    jac1$betadiver<-rep(0, nrow(jac1)) #create a column where the values of dissimilarity will be entered
    jac1<-jac1 %>% #with the function create at the beginning, create variables with the exact names of the pairwise comparison
      mutate(site1 = renamer(row, seq(1:nrow(.dados)), row.names(.dados)), 
             site2 = renamer(col, seq(1:nrow(.dados)), row.names(.dados))) %>%
      select(-row, -col) #remove the row and col columns
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
    jac1$betadiver <- as.vector(as.dist(jac0)) #put the dissimilarity values into the matrix of the run
    jac2 <- rbind(jac2, jac1) #put these results in the main matrix         
  }
  jac2[-1,] #before showing the matrix, remove the first line (as it was used just to create the base dataframe)
}