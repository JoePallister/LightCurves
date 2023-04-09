# Here we give code to generate some of the results from the PhD thesis 
# "Statistical inference on the periodicity of irregularly sampled light curves" by Efthymia Derezea

# We look at the generalized F-test and CVF test for hypothesis testing for periodic light curves, estimating the power of these tests. 

# We generate a number of periodic light curves using a GPR model and for each we generate a number of noises. Hence for each light curve 
# we construct a range of curves: light curve + noise. We call the ratio of the variance of the light curve to the variance of the noise 
# the signal to noise ratio (SNR). We fit a GPR model for these curves. The null hypothesis H_0 is that the data is just noise and the 
# alternative hypothesis H_1 is that the light curve has the chosen period (here we give the model M_1 the correct period). 

# We test the hypotheses using two different test statistics: generalized F and CVF. Both statistics' distributions under 
# the null hypothesis are not known, so we approximate them using saddlepoint approximation

# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------

#### Functions ####

require(MASS)

# The periodicic kernel function eq. (2.15)

kernel = function(t_i, t_j, A=1, h=1){
  sin_part = (-2/h^2)*(sin(pi*(t_i-t_j)/period))^2
  return (A*exp(sin_part))
}

# Matrix formed the kernel functions
kernel_mat = function(x_1, x_2, A=1, h=1){
  mat = matrix(1, nrow=length(x_1), ncol=length(x_2))
  for (i in 1:length(x_1)){
    for (j in 1:length(x_2)){
      mat[i, j] = kernel(t_i=x_1[i], t_j=x_2[j])
    }
  }
  return (mat)
}

# The function K and its derivatives, and the inverse of K' at 0. Page 47.
K = function(s, evalues){
  logs = log(1-2*s*evalues)
  return (-(1/2)*sum(logs))
}

K_prime = function(s){
  summands = evalues/(1-2*s*evalues)
  return (sum(summands))
}

K_prime_prime = function(s, evalues){
  summands = 2*evalues^2/(1-2*s*evalues)^2
  return (sum(summands))
}

K_prime_inverse_at_zero = function(){
  epsilon = 0.0000001
  lower_val = 1/(2*min_eig)
  upper_val = 1/(2*max_eig)
  # Sometimes, if all eigenvalues are negative, we have lower > upper so we need to reverse them
  if (lower_val > upper_val){
    temp = lower_val
    lower_val = upper_val
    upper_val = temp
  }
  root = uniroot(K_prime, lower=lower_val+epsilon, upper=upper_val-epsilon, tol = 1E-6)$root
  return (root)
}

# The saddlepoint approximation function. Page 46
F_saddlepoint= function(f_obs, M_1, evalues){
  s_hat = K_prime_inverse_at_zero()
  w = sign(s_hat)*(2*sqrt(abs(s_hat-K(s_hat, evalues))))
  u = s_hat*sqrt(abs(K_prime_prime(s_hat, evalues)))
  return (1-pnorm(w+(1/w)*log(u/w)))
}

# Here y, x_1 and x_2 are set later. These are not taken as parameters as we will want to optimize this function, but not over y, x_1 or x_2.
# Page 42
cross_valid_error = function(params){
  Ker = kernel_mat(x_1, x_2, h=param[1], A=param[3])
  temp = solve(Ker+params[2]^2*diag(length(x_1)))
  temp_diags = diag(temp)
  B_2 = diag(temp_diags)
  B = temp %*% B_2
  return (t(y) %*% B %*% t(B) %*% y)
}

#### end ####

# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------

#### Generalised F-test (sect 3.4.1) ####

# The period of our generated curves
period = 1.7

# The set of SNR values to loop over
SNR_set = seq(0.1, 1, len=20)

# The number of curves we generate for each SNR value
num_curves = 5

# The number of time points we use to draw curves
number_times = 200

# Select times from a uniform distribution
times = sort(runif(number_times, 1, 10))

# Empty matrix to collect results from loop
saddlepoint_results_F = matrix(0, ncol = num_curves, nrow = length(SNR_set))

# This kernel matrix gives the covariance for the multivariate normal sampling (2.15)
K_tt = kernel_mat(times, times)

# The loop
for (i in 1:length(SNR_set)){
  print(i)
  for (k in 1:num_curves){
    # Select number_times points from a multivariate normal with mean 0 and covariance given by the kernel matrix
    # The periodic form of the kernel ensures that times separated by a whole period give equal sampled values
    signal = mvrnorm(1, mu=rep(0, number_times), Sigma=K_tt)
    
    # Plot each sample
    # plot(times, signal)
    
    # Set seed so that noise changes with k
    set.seed(76454545+k)
    noise = rnorm(number_times, mean=1, sd=1)
    
    # Scale noise so that we have the correct SNR
    scale_factor = sqrt(var(signal))/(SNR_set[i]*var(noise))
    noise = scale_factor*noise
    
    # Full signal and errors
    full_signal = signal + noise
    e = full_signal - mean(full_signal)
    
    y = full_signal
    x_1 = times
    x_2 = times
    
    # Optimize the parameters in the cross-validation error (3.8)
    optim_params = optim(par =c(2.5,0.1,1), lower = c(0.1,0.07,0.1), upper=c(30,7.38,13), fn = cross_valid_error, method = "L-BFGS-B")$par
    optim_h = optim_params[1]
    optim_sigma = optim_params[2]
    optim_A = optim_params[3]
    
    # Matrices for calculating the p-values for the F-statistic (3.2). Page 40
    M_0 = diag(1, number_times) -(1/number_times) * matrix(1, ncol=number_times, nrow=number_times)
    # Page 25
    W = K_tt%*%solve(K_tt+optim_sigma^2*diag(1, number_times))
    # Page 40 
    M_1 = t(diag(1, number_times)-W)%*%(diag(1, number_times) - W)
    
    # We want to use (3.3) by replacing the y's in (3.2) but the transformation e = M_0y is not necessarily invertible. 
    # Apparently the following is a good approximation
    f_obs = (t(e)%*%(diag(1, number_times)-M_1)%*%e) / (t(e)%*%M_1%*%e)
    
    # We take eigenvalues (page 45) as follows
    evalues_mat = M_0%*%(diag(number_times)-(1+f_obs[1,1])*M_1)
    # This ensures that the quadratic form is represented by a symmetric matrix (hence real eigenvalues)
    evalues = eigen((1/2)*(evalues_mat + t(evalues_mat)), only.values=TRUE)$values
    min_eig = min(evalues)
    max_eig = max(evalues)
    # Approximate the p-value for the F-statstic
    saddlepoint_results_F[i, k] = F_saddlepoint(f_obs = f_obs, M_1 = M_1, evalues = evalues)
  }
}

sig_level = 0.02

# Probability of a curve having low enough p-value, for each SNR value. 
probs = rowSums(saddlepoint_results_F < sig_level) / num_curves
plot(SNR_set, probs, type='l')

#### end ####

# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------

#### CVF-test (sect 3.4.2) ####

# Note that almost all of this calculation is the same as for the generalized F-statistic, except for the  construction of the matrices M_0 and M_1

# The period of our generated curves
period = 1.7

# The set of SNR values to loop over
SNR_set = seq(0.1, 1, len=20)

# The number of curves we generate for each SNR value
num_curves = 5

# The number of time points we use to draw curves
number_times = 200

# Select times from a uniform distribution
times = sort(runif(number_times, 1, 10))

# Empty matrix to collect results from loop
saddlepoint_results_CVF = matrix(0, ncol = num_curves, nrow = length(SNR_set))

# This kernel matrix gives the covariance for the multivariate normal sampling (2.15)
K_tt = kernel_mat(times, times)

# The loop
for (i in 1:length(SNR_set)){
  print(i)
  for (k in 1:num_curves){
    # Select number_times points from a multivariate normal with mean 0 and covariance given by the kernel matrix
    # The periodic form of the kernel ensures that times separated by a whole period give equal sampled values
    signal = mvrnorm(1, mu=rep(0, number_times), Sigma=K_tt)
    
    # Plot each sample
    # plot(times, signal)
    
    # Set seed so that noise changes with k
    set.seed(76454545+k)
    noise = rnorm(number_times, mean=1, sd=1)
    
    # Scale noise so that we have the correct SNR
    scale_factor = sqrt(var(signal))/(SNR_set[i]*var(noise))
    noise = scale_factor*noise
    
    # Full signal and errors
    full_signal = signal + noise
    e = full_signal - mean(full_signal)
    
    y = full_signal
    x_1 = times
    x_2 = times
    
    # Optimize the parameters in the cross-validation error (3.8)
    optim_params = optim(par =c(2.5,0.1,1), lower = c(0.1,0.07,0.1), upper=c(30,7.38,13), fn = cross_valid_error, method = "L-BFGS-B")$par
    optim_h = optim_params[1]
    optim_sigma = optim_params[2]
    optim_A = optim_params[3]
    
    # Matrices for calculating the p-values for the CVF-statistic (3.9).
    # M_0 (3.10)
    temp = matrix(-1/(number_times-1), nrow=number_times, ncol=number_times)
    temp_2 = diag(1+1/(number_times-1), number_times)
    M = temp + temp_2
    M_0 = M^2
    
    # M_1, Page 42
    temp = K_tt+optim_sigma*diag(1, number_times)
    temp_2 = solve(temp)
    B_2 = diag(diag(temp_2))
    B = temp_2 %*% B_2
    M_1 = B %*% t(B)
    
    # We want to use (3.3) by replacing the y's in (3.2) but the transformation e = M_0y is not necessarily invertible. 
    # Apparently the following is a good approximation
    f_obs = (t(e)%*%(diag(1, number_times)-M_1)%*%e) / (t(e)%*%M_1%*%e)
    
    # We take eigenvalues (page 45) as follows
    evalues_mat = M_0%*%(diag(number_times)-(1+f_obs[1,1])*M_1)
    # This ensures that the quadratic form is represented by a symmetric matrix (hence real eigenvalues)
    evalues = eigen((1/2)*(evalues_mat + t(evalues_mat)), only.values=TRUE)$values
    min_eig = min(evalues)
    max_eig = max(evalues)
    # Approximate the p-value for the F-statstic
    saddlepoint_results_CVF[i, k] = F_saddlepoint(f_obs = f_obs, M_1 = M_1, evalues = evalues)
  }
}

# As before the following is the significance level
sig_level = 0.02

# Probability of a curve having low enough p-value, for each SNR value. 
probs = rowSums(saddlepoint_results_CVF < sig_level) / num_curves
plot(SNR_set, probs, type='l')

#### end ####