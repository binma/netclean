
silence<-function(mat){
  imat<-diag(70)
  dmat<-diag(diag((mat-imat)*mat))
  smat<-(mat-imat+dmat)*(1/mat)
  smat
}

data(mite)

mat<-as.matrix(1-vegdist(mite,diag=TRUE,upper=TRUE))