%% Questions 6 - Null Space and Range Space

% Load Data

load('mtxExamples.mat')

% Lets decompose them and then map to the null space...

[u1,s1,v1] = svd(mtx1) ;
[u2,s2,v2] = svd(mtx2) ;
[u3,s3,v3] = svd(mtx3) ;
[u4,s4,v4] = svd(mtx4) ;

s1_diag = diag(s1) ;
s2_diag = diag(s2) ;
s3_diag = diag(s3) ;
s4_diag = diag(s4) ;

mat1_null_index = (s1_diag < 10e-6) ;
mat2_null_index = (s2_diag < 10e-6) ;
mat3_null_index = (s3_diag < 10e-6) ;
mat4_null_index = (s4_diag < 10e-6) ;

mat1_range_index = (s1_diag > 10e-6) ;
mat2_range_index = (s2_diag > 10e-6) ;
mat3_range_index = (s3_diag > 10e-6) ;
mat4_range_index = (s4_diag > 10e-6) ;

% Matrix 2 and matrix 4 have null spaces.. lets transpose the input matrix
% and take a look at them.

MAT2_V_transposed = v2' ;
MAT4_V_transposed = v4' ;

mat_2_null = MAT2_V_transposed(mat2_null_index, :) ;
mat_4_null = MAT4_V_transposed(mat4_null_index, :) ;

% Now we can create random vectors which map to a zeros matrix

% For the second matrix - there are two vectors in the null space of the
% input matrix. Either of these vectors will work to give us a matrix that
% maps to the null space. We utilize this matrix, multiplied by  random
% scalar, to give us a random vector that maps to this zeros matrix.

mat2_null_rows  = mat_2_null((1:2),:) ;
random_mat2_vec_null = sum(mat2_null_rows * randn) ;

% Lets do this again for the fourth matrix...
% This time there is no need to utilize two vectors b/c there is only one
% in the nullspace

mat4_null_row  = mat_4_null(1,:) ;
random_mat4_vec_null = mat4_null_row * randn ;

% Now we can transform them with U,S,V^transposed

mat2_zeros_mat = u2*(s2*(v2'*random_mat2_vec_null'))
mat4_zeros_mat = u4*(s4*(v4'*random_mat4_vec_null'))

sprintf('All values are zero')

find(mat2_zeros_mat < 10e-6) 
find(mat4_zeros_mat < 10e-6) 

%% Mapping to the range space

% Take our Range Space vectors
mat_2_range = u2(:, mat2_range_index) ;
mat_4_range = u4(:, mat4_range_index) ;

% First we need to find the psuedo-inverse of S for each matrix.
% Take the transpose and then s diagonal values over themselves..
S_2_inv = s2' ;
S_2_inv(1,1) = 1/s2(1,1) ;    % There is only 1

% Then we need to isolate the input vector with this equation...
% ---> VS(pseudo)U'y(vec) = x(vec)

x_input_vec = v2*(S_2_inv*(u2'*(randn*mat_2_range))) ;
mtx2*x_input_vec 

% We can do the same for the second matrix, but it doesn't seem to work and
% it's 11:30..

S_4_inv = s4' ;

% First we need to find the psuedo-inverse of S for each matrix.
% Take the transpose and then s diagonal values over themselves..

S_4_inv(1,1) = 1/s4(1,1)  ; % There is only 1

x_input_vec2 = v4*(S_4_inv*(u4'*(randn*mat_4_range))) ;

mtx4*x_input_vec2

