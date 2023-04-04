
# Here we give code to generate some of the results from the PhD thesis 
# "Statistical inference on the periodicity of irregularly sampled light curves" by Efthymia Derezea

# In particular we are interested in the application of the standard F-test and generalised F-test to hypothesis testing
# for periodic light curves, so we estimate the power of these tests. Here the null hypothesis H_0 is that the data is just
# noise and the alternative hypothesis H_1 is that the light curve has the chosen period (here we give the model M_1 the correct period).  
# We generate a number of periodic light curves of the form a*sin + b*cos and for each we generate a number of noises. 
# Hence for each light curve we construct a range of curves: light curve + noise. We call the ratio of the variance of the light curve 
# to the variance of the noise the signal to noise ratio (SNR).

# For each light curve + noise we use two weighted linear models: y = const and y = a + b*sin + c*cos. H_0 is that the first model holds,
# meaning the function is just noise. H_1 is the second model, that we have a periodic function. We reject H_0 according to each
# of our F tests.

# For the generalised F-statistic the distribution under the null hypothesis is not known, so we approximate it using saddlepoint approximation. For this we need
# to invert a certain function K' (defined below). To see the behaviour of K' we give a detained example and perform an example saddlepoint
# approximation.

# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------

#### Standard F-statistic for light curve data ####

# The function we use to generate data. We will add this to some noise which will then be our simulated light curve
periodic_func = function(t, period=1, a1=1, a2=0) {
  return(a1*sin(2*pi*t/period) + a2*cos(2*pi*t/period))
}

# We construct a number of curves from a number of data points for each. We will loop over these curves
num_curves = 100
data_points = 200

# When constructing our light curves we will take different levels of noise, measured by the signal-noise ratio (SNR). For each curve we
# will loop over the SNR_set
SNR_set = seq(from=0.1, to=6, len=20)

# Initialize results matrix
results = matrix(0, nrow=length(SNR_set), ncol=1+num_curves)

# Set results col 1, the SNR values
results[,1] = SNR_set

# Select random sample points
times = 3*runif(data_points)

# Select random weights
weights = runif(data_points)

# Standardize weights
weights_standard = weights/sqrt(var(weights))

# The loop
for (i in 1:num_curves){
  if (i %% 25== 0){
    print('Curve')
    print(i)
  }
  # Take different period 1 curves, depending on parameters a1 and a2
  magnitude = periodic_func(times, period=1, a1=(i-1)/num_curves, a2=1-(i-1)/num_curves)
  for (j in 1:length(SNR_set)){
    SNR = SNR_set[j]
    # Scale weights based on SNR
    weights_scaled = weights_standard*sqrt(var(magnitude)/SNR)
    # Sample random errors
    errors = rnorm(data_points, sd=weights_scaled)
    # Add errors to generated data
    magnitude_and_errors = magnitude + errors
    # Linear models
    linreg_0 = lm(magnitude_and_errors ~ 1, weights = weights_scaled)
    linreg_1 = lm(magnitude_and_errors ~ sin(2*pi*times) + cos(2*pi*times), weights = weights_scaled)
    # Data to calculate F
    X_0 = model.matrix(linreg_0)
    X_1 = model.matrix(linreg_1)
    W = diag(weights_scaled)
    M_0 = diag(data_points) - X_0 %*% solve(t(X_0) %*% W %*% X_0) %*% t(X_0) %*% W
    M_1 = diag(data_points) - X_1 %*% solve(t(X_1) %*% W %*% X_1) %*% t(X_1) %*% W
    RSS_0 = t(magnitude_and_errors) %*% M_0 %*% magnitude_and_errors
    RSS_1 = t(magnitude_and_errors) %*% M_1 %*% magnitude_and_errors
    # Add F statistic to results
    results[j, 1+i] = RSS_0/RSS_1 - 1
  } 
}

results_frame = data.frame(results[,-1])

# We obtain the standard F-statistic by rescaling as follows
results_frame_standard_F = results_frame*(data_points-3)/(2)

# We reject the null hypothesis at a significance level=0.05, if values are larger than qf(0.95)
truth_standard_F = results_frame_standard_F > qf(0.95,data_points-3,2)

# Find probabilities of rejecting H_0
probs_standard_F = rowSums(truth_standard_F) / num_curves

# This gives the blue curve from the thesis (fig 3.4)
plot(SNR_set, probs_standard_F, xlab="SNR", ylab="Power")

#### end ####

# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------

#### Example of calculating the generalised F-statistic ####

# We calculate an example generalised F-statistic eq (3.2) using saddlepoint approximation, and we plot K' (which we will need to invert)
# to show the general behaviour

# Select some eigenvalues (in the full example these depend on the models M_0 and M_1 and the observed value of the F-statistic)
lambda_1 = -3
lambda_2 = 7
lambda_3 = 2

# We will need to invert K' (below) which has discontinuities at 1/(2*lambda), here these are -1/6, 1/14 and 1/4 
# The strange labelling of lambda gives a good labelling of mu below

# Set mu_i as the discontinuous values
mu_1 = 1/(2*lambda_1)
mu_2 = 1/(2*lambda_2)
mu_3 = 1/(2*lambda_3)

# The function K that we are interested in, and its first and second derivatives
K = function(s){
  return ((-1/2)*(log(1-2*s*lambda_1)+log(1-2*s*lambda_3)+log(1-2*s*lambda_3)))
}

K_prime = function(s){
  return (lambda_1/(1-2*s*lambda_1)+lambda_2/(1-2*s*lambda_2)+lambda_3/(1-2*s*lambda_3))
}

K_prime_prime = function(s){
  return (2*lambda_1^2/(1-2*s*lambda_1)^2+2*lambda_2^2/(1-2*s*lambda_2)^2+2*lambda_3^2/(1-2*s*lambda_3)^2)
}

# For a plot to see what is happening with K' we need to avoid plotting too close to the mu_i so we take some small parameters 
# which we will use to stay away from the mu_i
epsilon=0.3
epsilon_2=0.01

# We take 4 sets of values, corresponding to the 4 regions of K' separated by the discontinuities
values = seq(from=mu_1-epsilon, to=mu_1-epsilon_2, len=1000)
values_2 = seq(from=mu_1+epsilon_2, to=mu_2-epsilon_2, len=500)
values_3 = seq(from=mu_2+epsilon_2, to=mu_3-epsilon_2, len=500)
values_4 = seq(from=mu_3+epsilon_2, to=mu_3+epsilon, len=1000)

# Join the four sets
total = c(values, values_2, values_3, values_4)

# Plot of K'
plot(total, K_prime(total), xlab=expression(s), ylab=expression(K_prime(s)))
abline(v=mu_1, lty=2)
abline(v=mu_2, lty=2)
abline(v=mu_3, lty=2)

# We look for a solution to K'(s) = 0 between the largest negative mu and the smallest positive, here (-1/6, 1/14). In this range K' is 1:1
K_prime_inverse = function(){
  return (uniroot(K_prime, c(mu_1+epsilon_2, mu_2-epsilon_2), extendInt="upX",tol = 1E-6)$root)
}

# By definition s_hat, w and u are given by
s_hat = K_prime_inverse()
w = sign(s_hat)*(2*sqrt(abs(s_hat-K(s_hat))))
u = s_hat*sqrt(abs(K_prime_prime(s_hat)))

# Finally the approximation of the probability P(F>f_obs) is
F_stat = 1 -  pnorm(w+(1/w)*log(u/w))

#### end ####

# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------

#### Functions for saddlepoint approximation ####

# Here we collect the functions we will need below for the application of the saddlepoint approximation to generated light curves

# The generalised F-statistic eq. (3.2), with epsilon = M_0*y (page 38), gives, in the notation of Section 3.6, 
# A = M_0^T * (M_0 - M_1) * M_0 and B = M_0^T * M_1 * M_0, we want the eigenvalues of A - f_obs*B (page 45) which is
# M_0^T * (M_0-(1+t_obs)*M_1) * M_0. This may have complex eigenvalues so we need to pause

# A quadratic form in epsilon, written with a matrix C: epsilon^t * C * epsilon can be rewritten as: epsilon^T * ((1/2)*(C+C^T)) * epsilon
# without affecting the form. By doing this we have that (1/2)*(C+C^T) is symmetric so has real eigenvalues. Applying this to our problem
# we instead look at the eigenvalues of (1/2)*((A-f_obs*B) + (A-f_obs*B)^T), which are real. This is collected in evalues_func
evalues_func = function(f_obs, M_1, M_0){
  mat = t(M_0)%*%(M_0 - (1+f_obs)*M_1)%*%M_0
  values = eigen((1/2)*mat + (1/2)*t(mat), only.values=TRUE)$values
  return(values)
}

# The function K and its derivatives, and the inverse of K' at 0
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

# The saddlepoint approximation function
F_saddlepoint= function(f_obs, M_1, evalues){
  s_hat = K_prime_inverse_at_zero()
  w = sign(s_hat)*(2*sqrt(abs(s_hat-K(s_hat, evalues))))
  u = s_hat*sqrt(abs(K_prime_prime(s_hat, evalues)))
  return (1-pnorm(w+(1/w)*log(u/w)))
}

#### end ####

# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------

#### Saddlepoint approximation of F-statistic for light curve data ####

# For the light curve data the procedure is the same as for the standard F-statistic case until we calculate the generalised F-stastic

# Loop over different curves, all with period 1
num_curves = 100
data_points = 200

# Loop over SNR values
SNR_set = seq(from=0.1, to=1.5, len=10)

# Initialize results matrix
results = matrix(0, nrow=length(SNR_set), ncol=1+num_curves)

# Set results col 1, the SNR values
results[,1] = SNR_set

# Select random sample points
times = 3*runif(data_points)

# Select random weights from
weights = runif(data_points)

# Standardize weights
weights_standard = weights/sqrt(var(weights))

# The loop
for (i in 1:num_curves){
  if (i %% 25== 0){
    print(i)
  }
  # Take different period 1 curves, depending on parameters a1 and a2
  magnitude = periodic_func(times, period=1, a1=(i-1)/num_curves, a2=1-(i-1)/num_curves)
  for (j in 1:length(SNR_set)){
    SNR = SNR_set[j]
    # Scale weights based on SNR
    weights_scaled = weights_standard*sqrt(var(magnitude)/SNR)
    # Sample random errors
    errors = rnorm(data_points, sd=weights_scaled)
    # Add errors to generated data
    magnitude_and_errors = magnitude + errors
    # Linear models
    linreg_0 = lm(magnitude_and_errors ~ 1, weights = weights_scaled)
    linreg_1 = lm(magnitude_and_errors ~ sin(2*pi*times) + cos(2*pi*times), weights = weights_scaled)
    # Data to calculate generalised F-statstic
    X_0 = model.matrix(linreg_0)
    X_1 = model.matrix(linreg_1)
    W = diag(weights_scaled)
    M_0 = diag(data_points) - X_0 %*% solve(t(X_0) %*% W %*% X_0) %*% t(X_0) %*% W
    M_1 = diag(data_points) - X_1 %*% solve(t(X_1) %*% W %*% X_1) %*% t(X_1) %*% W
    RSS_0 = t(magnitude_and_errors) %*% M_0 %*% magnitude_and_errors
    RSS_1 = t(magnitude_and_errors) %*% M_1 %*% magnitude_and_errors
    # Here we diverge from the standard F code
    f_obs = (RSS_0/RSS_1 - 1)[1,1]
    # Get the eigenvalues we are interested in (see comments above the code for evalues_func)
    evalues = evalues_func(f_obs, M_1, M_0)
    min_eig = min(evalues)
    max_eig = max(evalues)
    results[j, 1+i] = F_saddlepoint(f_obs, M_1, evalues)
  } 
}

# Collect results
results_frame = data.frame(results[,-1])

# Reject H_0 if P < 0.002
reject = results_frame < 0.002 

# Find probabilities of rejecting H_0
probs_reject = rowSums(reject) / num_curves

# This gives the red curve from the thesis (figure 3.4)
plot(SNR_set, probs_reject)

#### end ####