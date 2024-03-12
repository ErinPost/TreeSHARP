// This Stan program defines a DTM regression model with shared horseshoe priors
// on the regression coefficients.
// This method is known as `treeSHARP`.

functions {

  real DM_lpmf( int [] y_v, vector gamma_v){
    int N = sum(y_v);
    real G = sum(gamma_v);
    return lgamma(G) - lgamma(N + G)
            + sum(lgamma(to_vector(y_v) + gamma_v))
            - sum(lgamma(gamma_v));
  }

  real DTM_lpmf(int [] y_i,
  vector gamma_i,
  int [] parents,
  int [] children,
  int [] n_children,
  int [] branches,
  int V){

    int pos = 1;
    vector[V] output;

    for(v in 1:V){
      int n = n_children[v];
      int C_v[n] = segment(children, pos, n);
      int B_v[n] = segment(branches, pos, n);
      int y_v[n] = y_i[C_v];
      vector[n] gamma_v = gamma_i[B_v];
      int N = sum(y_v);
      real G = sum(gamma_v);
    output[v] = lgamma(N + 1) + lgamma(G) - lgamma(N + G)
            + sum(lgamma(to_vector(y_v) + gamma_v))
           - sum(lgamma(gamma_v)) - sum(lgamma(to_vector(y_v) + 1));

    pos = pos + n;
    }
    return sum(output);
  }

}
// The input data is a vector 'y' of length 'N'.
data {
  int <lower=2> B; // number of branches for K category tree model
  int<lower=0> I; // number of groups  = nrow(X) = nrow(Y)
  int <lower=0>Y[I,B+1]; // observed breakdowns by branch
  int<lower=0> P ;//  number of reg parameters = ncol(X)
  matrix[I,P] X; // covariate matrix
  int V; // number of internal "parent" nodes
  int parents[V]; // list of parent nodes
  int children[B]; // list of child nodes
  int branches[B]; // list of branches that correspond to children above
  int n_children[V]; // number of children per parent
}

parameters {
  row_vector[B] alpha_raw ; // intercepts for each branch
  matrix[P,B] beta_raw;// regression paramters for each outcome
  vector<lower=0, upper = pi()/2 >[P] delta_unif ; // half cauchy instead of cauchy
  real<lower=0, upper = pi()/2 > phi_unif ;
}

transformed parameters {
  matrix[I,B] Alpha ;
  row_vector[B] alpha ;
  matrix[P,B] beta;
  matrix<lower=0>[I,B] gamma ;
  vector<lower=0>[P] delta;
  real<lower=0> phi ;

    phi = tan(phi_unif); // produces C+(0,1)... read as 0 + 1*tan(phi_unif)
    delta = tan(delta_unif);
    alpha = 10 * alpha_raw;
   for(i in 1:I){
     Alpha[i,] = alpha ; // 10 * alpha because SD = 10
   }

   for(p in 1:P){
     for(b in 1:B){
       beta[p,b] = phi * delta[p] * beta_raw[p,b];//delta[p]normal(0,delta[p] * phi)
     }
   }

  gamma = exp(Alpha + X * beta);


}

model {
for(p in 1:P){
  for(b in 1:B){
    beta_raw[p,b] ~ normal(0,1) ;
  }
}

alpha_raw ~ normal(0,1);
phi_unif ~ uniform(0,pi()/2);
delta_unif ~ uniform(0,pi()/2);

for(i in 1:I){
  vector[B] gamma_i = to_vector(gamma[i,]);
  target += DTM_lpmf(Y[i,]|gamma_i, parents, children,n_children,branches,V);
}

}

generated quantities {

    vector[I] log_lik;
  for (i in 1:I) {
    vector[B] gamma_i = to_vector(gamma[i,]);
    log_lik[i] = DTM_lpmf(Y[i,]|gamma_i, parents, children,n_children,branches,V);
  }

}






