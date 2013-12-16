decon<-function(mat){
beta = 0.6;
alpha = 0.8;
control=0;

n = dim(mat)[1];
mat = mat*(1-diag(n));

#% thresholding the input matrix

y=quantile(c(mat),1-alpha);
mat_th=mat*(mat[mat>=y]);
       
         #***********************************
          #% eigen decomposition
          #%disp('decomposition and deconvolution...')

U<-eigen(mat_th)$vectors
D<-diag(eigen(mat_th)$values)         
          
lam_n=abs(min(min(diag(D)),0));
lam_p=abs(max(max(diag(D)),0));
          
m1=lam_p*(1-beta)/beta;
m2=lam_n*(1+beta)/beta;
m=max(m1,m2);
          
          
         # %network deconvolution
for(i in 1:dim(D)[1]){
          D[i,i] = (D[i,i])/(m+D[i,i]);
}
mat_new1 = U*D*t(U);

ind_edges = (mat_th[mat_th>0])*1.0;
ind_nonedges = (mat_th[mat_th==0])*1.0;
m1 = max(max(mat*ind_nonedges));
m2 = min(min(mat_new1));
mat_new2 = (mat_new1+max(m1-m2,0))*ind_edges+(mat*ind_nonedges);
m1 = min(min(mat_new2));
m2 = max(max(mat_new2));
mat_nd = (mat_new2-m1)/(m2-m1);

}