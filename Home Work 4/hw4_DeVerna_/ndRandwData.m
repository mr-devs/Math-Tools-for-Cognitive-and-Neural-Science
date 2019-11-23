function outData = ndRandwData(data, mean, cov)
% This function takes in a dataset of unit vectors and a desired mean and
% covariance matrix with which to transform that data by.
%   'data' = the unit vecot data to transform
%   'mean' = N-vector
%   'cov'ariance = NxN matrix
% OUTPUT = Your 'data' transformed from its original shape 
% to match the mean and covariance that you input to the function. 
% Dimensiosn = num rows & N columns

[V, D] = eig(cov)                   ; % We take the eigen vectors (V) and the eigenvalues (D) to apply to M (eventually)
M = V * sqrt(D)                     ; % Create the core of the trasformation matrix (gives you the scale and rotation needed)
outData = M*data + mean'            ; % Create fully transformed data
outData = outData'                  ;

end

