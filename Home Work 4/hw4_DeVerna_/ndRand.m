function samples = ndRand(mean, cov, num)
% Generates a set of samples drawn from an N-dimensional Gaussian
% distribution with the specified parameters:
%   'mean' = N-vector
%   'cov'ariance = NxN matrix
%   'num' = [optional -> default = 1] specifies # of samples to return
% OUTPUT = A sample of data points transformed from a Gaussian Distribution
% to match the mean and covariance that you input to the function. 
% Dimensiosn = num rows & N columns

data = randn([length(mean) num])    ;

[V, D] = eig(cov)                   ; % We take the eigen vectors (V) and the eigenvalues (D) to apply to M (eventually)
M = V * sqrt(D)                     ; % Create the core of the trasformation matrix (gives you the scale and rotation needed)
samples = M*data + mean'            ; % Create fully transformed data, need to add mean b/c it's currently centered around zero
samples = samples'                  ;

end

