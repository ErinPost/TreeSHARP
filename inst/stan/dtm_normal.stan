
functions {

  // DTM likelihood is a product of DM likelihoods for each parent node
  // stackoverflow.com/questions/58191561/matrix-with-simplex-columns-in-stan
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
    //int V = num_elements(parents);
    vector[V] output;
    //vector C_v;
    for(v in 1:V){
      int n = n_children[v];
      int C_v[n] = segment(children, pos, n);
      int B_v[n] = segment(branches, pos, n);
      int y_v[n] = y_i[C_v];
      vector[n] gamma_v = gamma_i[B_v]; //segment(gamma_i, pos, n);
      int N = sum(y_v);
      real G = sum(gamma_v);
    output[v] = lgamma(N + 1) + lgamma(G) - lgamma(N + G)
            + sum(lgamma(to_vector(y_v) + gamma_v))
           - sum(lgamma(gamma_v)) - sum(lgamma(to_vector(y_v) + 1));
  //output[v] = -0.0888;

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
  int branches[B];
  int n_children[V]; // number of children per parent
  real <lower = 0> sigma ; // SD of normal dist for betas
}


parameters {
  row_vector[B] alpha ; // intercepts for each branch
  matrix[P,B] beta ;// regression paramters for each outcome
}

transformed parameters {
  matrix[I,B] Alpha ;
  matrix<lower=0>[I,B] gamma ;
   for(i in 1:I){
     Alpha[i,] = alpha ;//#exp(alpha + X[i,] * beta) ;
   }
  gamma = exp(Alpha + X * beta);
}

model {


for(p in 1:P){
beta[p,] ~ normal(0,sigma); // switch back to 3
}
alpha ~ normal(0, 10);
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







