function convolvedMat = createConvMat(x,M)
%createConvMat Creates a convolution matrix based on an input vector and
%given dimensionality. The resulting matrix 'convolvedMat' allows you to do convolution
%between the input vector 'x' and some other vector 'v' by multiplying 'v'
%by the 'convolvedMat' matrix and not using the conv() function.

% INPUTS:
% x = input vector
% M = some dimensionality

% Find the important dimension
N = length(x) + M - 1;

% Make a new matrix with the correct final matrix dimensions. Each is a
% copy of the original matrix with zeros afterwards
copies_w_zeros = repmat([x(:); zeros(M,1)], 1, M+1);

% Now we turn this matrix into a very long row vector, and reshape it into the
% proper dimensions. This gives us the result we want.
convolvedMat = reshape(copies_w_zeros(1:N*M), N, M);

% The above returned answer would be incorrect if a row vector is
% originally input(because it forces 'dimension' to equal rows). 
% We can correct this by transposing the final result if x is a row.
if isrow(x)
    convolvedMat = convolvedMat';
end

end

