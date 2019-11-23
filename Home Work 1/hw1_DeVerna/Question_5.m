%% Gram-Schmidt

%{
A classic method for constructing an orthonormal basis is known as Gram-
Schmidt orthogonalization. 
- First, one generates an arbitrary unit vector (e.g., by normalizing a
vector created with randn). 
- Each subsequent basis vector is created by generating another arbitrary vector, 
subtracting off the projections of that vector along each of the previously
created basis vectors, and normalizing the remaining vector.

- Write a MATLAB function gramSchmidt that takes a single argument, N, specifying the
dimensionality of the basis. It should then generate an N×N matrix whose columns contain
a set of orthogonal normalized unit vectors. 
- Try your function for N = 3, and plot the basis vectors (you can use MATLAB’s 
rotate3d to interactively examine these). Check your function numerically by 
calling it for an N larger than 3 and verifying that the resulting matrix is 
orthonormal. 
%}

% I could not figure out the recursive portion!

n = 3 ;

gramSchmidt(n)


