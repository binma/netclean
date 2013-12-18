deconreg<-function(mat,beta=0.95,alpha=1,control=0){}

#**********************
# removing self-loops

n_tf=dim(mat)[1];
for (i in 1:n_tf){
mat(i,i)=0;
}

#**********************
# making the TF-TF network symmetric
# note some algorithms only output one-directional edges
# the other direction is added in that case

tf_net=mat[1:n_tf,1:n_tf];
ptf=which(tf_net!=t(tf_net));
   xx=ptf[,1];yy=ptf[,2];          
             tf_net_final=tf_net;
             for (i in 1:length(xx)){
             
             if (tf_net[xx[i],yy[i]]!=0 & tf_net[yy[i],xx[i]]!=0)
             tf_net_final[xx[i],yy[i]]= [tf_net[xx[i],yy[i]]+tf_net[yy[i],xx[i]]]/2;
             tf_net_final[yy[i],xx[i]]=tf_net_final[xx[i],yy[i]];
             else if (tf_net[xx[i],yy[i]]==0)
             tf_net_final[xx[i],yy[i]]= tf_net[yy[i],xx[i]];
             tf_net_final[yy[i],xx[i]]=tf_net_final[xx[i],yy[i]];
             else if(tf_net[yy[i],xx[i]]==0)
             tf_net_final[xx[i],yy[i]]= tf_net[xx[i],yy[i]];
             tf_net_final[yy[i],xx[i]]=tf_net_final[xx[i],yy[i]];
           
             }
             
             mat[1:n_tf,1:n_tf]=tf_net_final;
             
             #**********************
             # setting network density to alpha
             
             y=quantile(c(mat),1-alpha);
             mat_th=mat*(mat[mat>=y]);
             
             mat_th[1:n_tf,1:n_tf]=(mat_th[1:n_tf,1:n_tf]+t(mat_th[1:n_tf,1:n_tf]))/2;
temp_net=(mat_th>0)*1.0;
temp_net_remain=(mat_th==0)*1.0;
mat_th_remain=mat%*%temp_net_remain;
m11=max(max(mat_th_remain));

#**********************
  # check if perturbation is needed

if (control_p!=1)
# padding zero to make it a square matrix
mat1=[mat_th;zeros(n-n_tf,n)];
mat1=full(mat1);
# decomposition step
[U,D] = eigen(mat1);
if rcond(U)<10^-10
control_p=1;
end
end
#**********************
  # perturb input slightly
if control_p==1
r_p=0.001;
rand_tf=r_p*rand(n_tf);
rand_tf=(rand_tf+rand_tf')/2;
         for i=1:n_tf
         rand_tf(i,i)=0;
         end
         
         rand_target=r_p*rand(n_tf,n-n_tf);
         mat_rand=[rand_tf,rand_target];
         mat_th=mat_th+mat_rand;
         
         #******
         # padding zero to make it a square matrix
         mat1=[mat_th;zeros(n-n_tf,n)];
         mat1=full(mat1);
         #******
         # decomposition step
         [U,D] = eig(mat1);
         end
         
         #******
         # scaling based on eigenvalues
         
         lam_n=abs(min(min(diag(D)),0));
         lam_p=abs(max(max(diag(D)),0));
         #
         m1=lam_p*(1-beta)/beta;
         m2=lam_n*(1+beta)/beta;
         scale_eigen=max(m1,m2);
         # applying network deconvolution filter
         for i=1:n
         D(i,i)=(D(i,i))/(scale_eigen+D(i,i));
         end
         
         net_new=U*D*inv(U);
         #**********************
         # adding remaining edges
         
         net_new2=net_new(1:n_tf,:);
         m2=min(min(net_new2));
         net_new3=(net_new2+max(m11-m2,0)).*temp_net;
         mat_nd=net_new3+mat_th_remain;
         
         
         