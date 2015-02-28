#this function rarefy the entire dataset to the lower species richness observed
#the user should provide this levels, so that the function knows which value is it
#although this is the easiest approach, you may lose important information between
#pairs of species rich communities, as they're the ones that will lose most of the
#species in the rarefying procedure


beta.rare<-function(.dados, .reps, .comp) {
  reps<-.reps
  #create a data frame to store the information
  jac0<-as.data.frame(matrix(0, ncol=nrow(.dados), nrow=nrow(.dados)))
  for(i in 1:reps) {
    #rarefy the entire dataset
    r1<-rrarefy(.dados, sample = .comp)
    #store the results of vegdist as a dataframe
    jac1<-as.data.frame(as.matrix(vegdist(r1, method="jac")))
    #and sum the results for each run
    jac0<-jac0+jac1
  }
  r2<-jac0/reps #take the average of that dissimilarity
  rownames(r2)<-rownames(.dados) #fix rownames
  colnames(r2)<-rownames(.dados) #fix colnames
  as.dist(r2) #show the output
}