functions {
  real fiduciaJacobian_lp(int n, int d, vector y, matrix X, vector beta, real sigma2a, real sigma2e, matrix Sa, matrix Se, vector weight){
    matrix[n,n] S;
    matrix[n,d+2] fiducialM;
    
    S=sigma2a*Sa+sigma2e*Se;
    fiducialM=diag_pre_multiply(weight, append_col(X,append_col(Se*mdivide_left_spd(S,y-X*beta),Sa*mdivide_left_spd(S,y-X*beta))));
    return 0.5*log_determinant(fiducialM'*fiducialM)-0.5*(d+2)*log(n);
  }
}

data {
  int<lower=1> n; //number of cases
  int<lower=1> d; //number of predictors
  real<lower=0> k; //max eigenvalue of Sa - typically max number of replication in a level
  
  matrix[n,d] X; //design matrix
  matrix[n,n]  Sa;  //random effect covariance matrix (typically block diagonal taking values 0,1)
  vector[n] y; //dependent variables
}

transformed data{
  cov_matrix[n]  Se;  //random effect covariance matrix - identity matrix
  vector[n] unitv;
  vector<lower=0>[n] weight;
  for(i in 1:n){
    unitv[i]=1;
    weight[i]=1/sqrt(sum(Sa[i]));
  }
  Se=diag_matrix(unitv);
}

// The parameters accepted by the model. Our model
parameters {
  vector[d] beta; //fixed effects
  real<lower=0> sigma2e; //residual variance
  real<lower=0> u2; //sigma2e + k*sigma2a = k*u2 >0 
}

transformed parameters{
  real sigma2a; //random effect variance
  
  sigma2a=u2-sigma2e/k;
}

model {
  y ~ multi_normal(X*beta, sigma2a*Sa+sigma2e*Se);
  target+=fiduciaJacobian_lp(n, d, y, X, beta, sigma2a, sigma2e, Sa, Se, weight);
}


