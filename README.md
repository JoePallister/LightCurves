# LightCurves
We are interested in methods of period detection for light curves

In "FTestsWeightedLinearRegression.R" we give code to generate some of the results from the PhD thesis "Statistical inference on the periodicity of irregularly 
sampled light curves" by Efthymia Derezea

In particular we are interested in the application of the standard F-test and generalised F-test to hypothesis testing for periodic light curves, so we estimate 
the power of these tests. Here the null hypothesis H_0 is that the data is just noise and the alternative hypothesis H_1 is that the light curve has the chosen 
period (here we give the model M_1 the correct period). We generate a number of periodic light curves of the form a*sin + b*cos and for each we generate a number 
of noises. Hence for each light curve we construct a range of curves: light curve + noise. We call the ratio of the variance of the light curve to the variance 
of the noise the signal to noise ratio (SNR).

For each light curve + noise we use two weighted linear models: y = const and y = a + b*sin + c*cos. H_0 is that the first model holds, meaning the function is just 
noise. H_1 is the second model, that we have a periodic function. We reject H_0 according to each of our F tests.

For the generalised F-statistic the distribution under the null hypothesis is not known, so we approximate it using saddlepoint approximation. For this we need
to invert a certain function K' (defined below). To see the behaviour of K' we give a detained example and perform an example saddlepoint approximation.
