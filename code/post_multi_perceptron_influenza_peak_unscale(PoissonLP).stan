functions {
 matrix calculate_alpha(matrix X, matrix beta_first, matrix[] beta_middle, matrix beta_last, matrix weight_first, matrix[] weight_middle, matrix weight_last, int N, int num_hidden_layers, int num_hidden_nodes){
    matrix[rows(X), cols(beta_first)] layer_values_first;
    matrix[rows(X), cols(beta_first)] layer_values_middle[num_hidden_layers];
    matrix[rows(X), cols(beta_last)] alpha;
    layer_values_first = tanh(beta_first + X * weight_first);   
    layer_values_middle[1] = tanh(beta_middle[1] + layer_values_first * weight_middle[1]);
    for(i in 2:(num_hidden_layers)){
      layer_values_middle[i] = tanh(beta_middle[i] + layer_values_middle[i-1] * weight_middle[i]);
    }
    alpha = beta_last + layer_values_middle[num_hidden_layers] * weight_last;
    return alpha;
  }

 matrix calculate_peak(matrix alpha_post_2){
   matrix[rows(alpha_post_2), cols(alpha_post_2)] alpha_post_2_trans;
   for(i in 1:cols(alpha_post_2)){
     if (alpha_post_2[1,i] <= 183)
      alpha_post_2_trans[1,i] = 50;
     else
      alpha_post_2_trans[1,i] = alpha_post_2[1,i];
   }
   return alpha_post_2_trans;
 }
 
 matrix unscale_alpha(matrix alpha_post, matrix y_maxs, matrix y_mins){
   matrix[rows(alpha_post), cols(alpha_post)] alpha_post_unscaled;
   matrix[rows(alpha_post), 1] vec_ones;
   matrix[rows(alpha_post), cols(alpha_post)] mat_mins;
   matrix[rows(alpha_post), cols(alpha_post)] mat_maxs;
   
   vec_ones = rep_matrix(1, rows(alpha_post), 1);
   mat_mins = vec_ones * y_mins;
   mat_maxs = vec_ones * y_maxs;
   alpha_post_unscaled = (alpha_post + 1) .* (mat_maxs - mat_mins) / 2 + mat_mins;
   
   return alpha_post_unscaled;
}

 real log_poisson_lpdf(vector y_targets_1, vector alpha_post_1){
   real lp = 0;
   for (i in 1:cols(y_targets_1)){
     lp += y_targets_1[i]*log(alpha_post_1[i]) - alpha_post_1[i];
   }
   return lp;
   }

}

data {
  int<lower=0> num_targets;
  int<lower=0> num_inputs;
  int<lower=0> num_outputs;
  int<lower=0> num_hidden_nodes;
  int<lower=1> num_hidden_layers;
  matrix[1, num_outputs] y_maxs;
  matrix[1, num_outputs] y_mins;
  matrix[num_targets,num_outputs] y_targets;
  //matrix[num_targets,num_outputs] y_targets_se;
  matrix[num_inputs, num_hidden_nodes] weight_first;
  matrix[num_hidden_nodes, num_hidden_nodes] weight_middle[num_hidden_layers];
  matrix[num_hidden_nodes, num_outputs] weight_last;
  matrix[1, num_hidden_nodes] beta_first;
  matrix[1, num_hidden_nodes] beta_middle[num_hidden_layers];
  matrix[1, num_outputs] beta_last;
}

transformed data{
  matrix[1, 16] y_targets_1;
  matrix[1, 4] y_targets_2;
//  matrix[1, 16] y_targets_se_1;
//  matrix[1, 4] y_targets_se_2;
  
  y_targets_1 = y_targets[:,1:16];
  y_targets_2 = y_targets[:,17:20];
  //y_targets_se_1 = to_matrix(y_targets_se[1,1:16]);
  //y_targets_se_2 = to_matrix(y_targets_se[1,17:20]);  
}

parameters {
  matrix<lower=-1, upper=1>[num_targets,num_inputs] Xq;
}


model{
  matrix[1, num_outputs] alpha_post;
  matrix[1, num_outputs] alpha_post_unscaled;
  matrix[1, num_outputs] alpha_post_unscaled_exp;
  matrix[1, 16] alpha_post_1;
  matrix[1, 4] alpha_post_2;
  matrix[1, 4] alpha_post_2_trans;
    alpha_post = calculate_alpha(Xq, beta_first, beta_middle, beta_last, weight_first, weight_middle, weight_last, num_targets, num_hidden_layers, num_hidden_nodes);
    // unscale the ANN results
    alpha_post_unscaled = unscale_alpha(alpha_post, y_maxs, y_mins);
    // exponentiate the unscaled ANN outcomes
    for (i in 1:num_outputs){
      alpha_post_unscaled_exp[1,i] = exp(alpha_post_unscaled[1,i]);
    }
    //break down 16 vs 4 outcomes
    alpha_post_1 = alpha_post_unscaled_exp[:,1:16];
    alpha_post_2 = alpha_post_unscaled_exp[:,17:20];
    alpha_post_2_trans = calculate_peak(alpha_post_2);
    //calculate Poisson deviance
    target += sum(to_vector(y_targets_1) .* log(to_vector(alpha_post_1)) - to_vector(alpha_post_1));
    //to_vector(y_targets_1) ~ lognormal(to_vector(alpha_post_1),to_vector(y_targets_se_1)); //get SE directly from data
    to_vector(y_targets_2) ~ normal(to_vector(alpha_post_2_trans),10); //get SE directly from data
    to_vector(Xq) ~ uniform(-1,1);
    //to_vector(Xq) ~ normal(0,1);
}




