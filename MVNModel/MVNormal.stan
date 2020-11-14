functions {
  matrix myIplusA(int d, vector aa){
    matrix[d,d] IplusA;
    int icount;
    
    icount=1;
    for(i in 1:d){
      for(j in 1:i){
        if(i==j) IplusA[i,j]=1;
        else{
          IplusA[i,j]=aa[icount];
          IplusA[j,i]=-aa[icount];
          icount=icount+1;
        }
      }
    }
    return IplusA;
  }
  
  real fiduciaJacobian_lp(int n, int d, matrix SSE, matrix IplusA, matrix U){
    matrix[d,d] IplusAinv;
    matrix[d,d] UtSSE;
    matrix[d,d] IplusAinvSSE;
    matrix[(d*(d+1))/2,(d*(d+1))/2] JacobianMatrix;
    //The derivatives with respect to mu do not contribute
    int icount; //counters for populating the covarianca matrix
    int jcount;
    
    UtSSE=U'*SSE;
    IplusAinv=inverse(IplusA);
    IplusAinvSSE=(IplusA')\SSE;
    
    icount=0;
    for(k in 1:d){    //These nested for loops are for iterating
      for(l in 1:k){  //through the rows
        icount+=1;
        jcount=0;
        for(i in 1:d){   //These nexted for loops are for iterating
          for(j in 1:i){//through the cols
            jcount+=1;
            if(i==j){//Derivative with respect to lambda
              JacobianMatrix[icount,jcount]=U[k,i]*UtSSE[i,l]+U[l,i]*UtSSE[i,k];
            }
            else{
              JacobianMatrix[icount,jcount]=-2*(IplusAinv[k,i]*IplusAinvSSE[j,l]-IplusAinv[k,j]*IplusAinvSSE[i,l]+IplusAinv[l,i]*IplusAinvSSE[j,k]-IplusAinv[l,j]*IplusAinvSSE[i,k]);
            }
          }
        }
      }
    }
    
    return log_determinant(JacobianMatrix) - 0.5*d*(d+1)*log(n);
  }
}

data {
  int<lower=1> n;
  int<lower=1> d;
  vector[d] ybar;
  cov_matrix[d] SSE;
}

parameters {
  vector<lower=-1,upper=1>[(d*(d-1))/2] aa;
}

transformed parameters{
  matrix[d,d] IplusA;
  matrix[d,d] U;
  vector[d] dUtSnU;
  real logJacA;
  real logMarginalLikelihood;
  
  IplusA=myIplusA(d,aa);
  U=IplusA'/IplusA;
  dUtSnU=diagonal(U'*SSE*U);
  
  logMarginalLikelihood=-0.5*(n-1)*sum(log(dUtSnU))+d*lgamma(0.5*(n-1))-d*log(2)-0.5*d*(n-1)*log(pi())-0.5*d*log(n);
  logJacA=fiduciaJacobian_lp(n, d, SSE, IplusA, U);
}

model {
  aa~uniform(-1,1); //I am not sure if this is needed
  
  target+=logMarginalLikelihood; //likelihood + the lamda part of Jacobian
  target+=logJacA; //adjusting for the rest of the Jacobian
}

generated quantities{
  cov_matrix[d] Sigma;
  cov_matrix[d] invSigma;
  vector<lower=0>[d] lambda; 
  vector[d] mu;
  
  for(i in 1:d){
    lambda[i]=sqrt(1/gamma_rng(0.5*(n-1),0.5*dUtSnU[i]));
  }
  
  Sigma=diag_post_multiply(U,lambda)*diag_post_multiply(U,lambda)';
  invSigma=diag_post_multiply(U,1.0./lambda)*diag_post_multiply(U,1.0./lambda)';
  
  mu=multi_normal_rng(ybar,Sigma/n);
}
