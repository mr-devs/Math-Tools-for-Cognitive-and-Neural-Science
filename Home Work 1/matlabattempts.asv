%% Math Tools I - HW #1

% Purpose: This code has been written to wrestle with the first homework
% problems presented in the infamous Math Tools course at New York
% University.

% Author: Matthew DeVerna

% Date: 09/24/19

%% Question 1: Determine if a system is linear. 
% If it is, create a matrix which satsifiys the input/output pairs.

% In order for a system to be linear it must meet the requirements of
% superposition: homogeneity and additivity. More simply, any type of linear
% combination to the inputs should have the same effect on the
% outputs.
    
%% Question 1: System 1:

input1 = 1          ;
input2 = 2.5        ;
output1 = [4;6]     ;
output2 = [10;14]   ;

% Because scaling input1 by 2.5 gives input2 exactly, the same affect
% should be seen for the outputs - if the system is linear.

if input1 * 2.5 == input2                           % Set conditional if-statement
    disp('Input 1 * 2.5 equals Input 2')            % If the condition is accurate, print this.

else
    disp('Input 1 * 2.5 DOES NOT equal Input 2')    % If the condition is NOT accurate, print this.
end

% Lets try the same for the outputs...


if output1 * 2.5 == output2                         % Set conditional if-statement
    disp('Output 1 * 2.5 equals Output 2')          % If the condition is accurate, print this.
else
    disp('Output 1 * 2.5 DOES NOT equal Output 2')  % If the condition is NOT accurate, print this.
end

% As we can see, this system does not follow the laws of superposition
% (specifically homogeneity). Thus, we conclude that this system is not linear.

%% Question 1: System 2:

% Since the scope of this questions programming complexity is quite simple,
% I will use the same variables as the last question. The linear nature of
% this script allows me to do so without worrying about workspace/variable
% issues - and offers the benefit of simple, clear variables.

input1 = [6;3]      ;
input2 = [-2;-1]    ;
output1 = [12;12]   ;
output2 = [-6;-6]   ;

% Because scaling input1 by -1/3 gives input2 exactly, the same affect
% should be seen for the outputs - if the system is linear.

if input1 * (-1/3) == input2           % Set conditional if-statement
    disp('Input 1 * -1/3 equals Input 2')            % If the condition is accurate, print this.

else
    disp('Input 1 * -1/3 DOES NOT equal Input 2')    % If the condition is NOT accurate, print this.
end

% Lets try the same for the outputs...

if output1 * (-1/3) == output2         % Set conditional if-statement
    disp('Output 1 * -1/3 equals Output 2')          % If the condition is accurate, print this.
else
    disp('Output 1 * -1/3 DOES NOT equal Output 2')  % If the condition is NOT accurate, print this.
end

% Once again, we can see that this system does not follow the laws of superposition
% (specifically homogeneity). Thus, we conclude that this system is not linear.

%% Question 1: System 3

input1 = [1;2]  ;
input2 = [1;-2] ;
input3 = [3;0]  ;
output1 = [5;-1];
output2 = [1;4] ;
output3 = [7;8] ;

%{ 
Because scaling input3 by 2/3, and then subtracting input 1 from that
product results in input2 exactly, the same affect should be seen for 
the outputs - if the system is linear.
%}

if (input3 * 2/3) - input1 == input2                        % Set conditional if-statement
    disp('Input3 * 2/3) - Input1 equals Input 2')           % If the condition is accurate, print this.

else
    disp('Input3 * 2/3) - Input1 DOES NOT equals Input 2')  % If the condition is NOT accurate, print this.
end

% Lets try the same for the outputs...

if (output3 * 2/3) - output1 == output2                         % Set conditional if-statement
    disp('Output 3 * 2/3 - Output1 equals Output 2')            % If the condition is accurate, print this.
else
    disp('Output 3 * 2/3 - Output1 DOES NOT equals Output 2')   % If the condition is NOT accurate, print this.
end

% Daaamn, my man, these systems just DO NOT want to follow the laws of superposition
% (specifically homogeneity). Thus, we conclude that this system is not linear.

%% Question 1: System 4

input1 = [2 ; 4];
input2 = [-2; 1];
output1 = 0 ;
output2 = 3 ;

%{ 
Based on the values above, you can tell immediately that this system
cannot be linear. Starting with the outputs, there is nothing that you can scale 
output1 (zero) by to get to output2 (three). Additionally, in order to get an output
of zero from input1, we can solve this problem by taking the dot product of
input1 with a 'testVec[tor]' of [1,(-1/2)]. If this is applied to input2, 
we find that this does not result in output2. 
%}

% Create our test vector...
testVec = [1,(-1/2)] ;

if (testVec * input1) == output1                        % Set conditional if-statement
    disp('The dot product of testVec and Input1 equals Output1')           % If the condition is accurate, print this.

else
    disp('The dot product of testVec and Input1 DOES NOT equals Output1')  % If the condition is NOT accurate, print this.
end

% Lets try the same for the input2 and output2...

if (testVec * input2) == output2                        % Set conditional if-statement
    disp('The dot product of testVec and Input2 equals Output1')           % If the condition is accurate, print this.

else
    disp('The dot product of testVec and Input2 DOES NOT equals Output1')  % If the condition is NOT accurate, print this.
end

% Additionally, if we want to try and scale input1 to input2, we would use
% the vector [-1,(1/4)].

[1,(1/4)]' * input1


% We should be able to add some factor to both the inputs and outputs, and
% get the same answer.

%% Question 1: System 5

input1 = 0;
ouput1 = [1;2];

% We can tell right away that this is not linear because we are starting at
% the origin. By definition, there is no linear combination that can take
% us away from the origin and add a new dimension. 


%% Question 2: Inner Product with a Unit Vector

%{
Given unit vector u and arbitrary vector v - write expressions for
computing:

--> PART 1
a) the component of v lying along the direction of u,
(b) the component of v that is orthogonal (perpendicular) to u, and
(c) the distance from v to the component that lies along direction u.

%}


% Create random vectors
v = randn(1,2) ;
u = randn(1,2) ;

% (a) calculate the component of v lying along the direction of u
lenVec2 = sqrt(u*u');
unitVec2 = u/lenVec2 ;
newLineMagnitude = v*unitVec2' ;
projectionUnit = newLineMagnitude * unitVec2 ;

% According to the lab notes, it should just be the residual...
residual = v - projectionUnit ;
orthogonalLine = projectionUnit + residual;

% Open a figure
projectionThingy = figure ;
hold on

% Create names for each line to plot (for the legend call below)
name1 = 'Original Vector'   ;  
name2 = 'Projection'        ;
name3 = 'Unit Vector'       ;
name4 = 'Residual'          ;

% Plot the original vector 
vPlot           = plot([0,v(1)],[0,v(2)],...
                   'b', 'LineWidth',2) ; 

% Plot the projection vector
projectionPlot  = plot([0,projectionUnit(1)],[0,projectionUnit(2)],...
                  'r', 'LineWidth',2); name3 = 'Projection';

% Plot the unit vector
unitPlot        = plot([0,unitVec2(1)],[0,unitVec2(2)],...
                    'k--*', 'LineWidth',2); name3 = 'Unit Vector';

% Plot the unit vector
residualPlot    = plot([projectionUnit(1), orthogonalLine(1)], [projectionUnit(2),orthogonalLine(2)],...
                    'g', 'LineWidth',2); name4 = 'Residual';
% Set axis to equal proportions
axis equal

legend(name1, name2, name3, name4, 'location', 'northeastoutside')
title('Projection/Orthogonal Component Plot')
xlabel('x-axis')
ylabel('y-axis')

%% Question 2: PART 2--> Proving this works for higher dimensions

%{
This section seeks to prove the above section in higher dimensions by
writing code to verify...
1) the vector in (a)- the projection line - is along the same line as the
unit vector
2) 
 
The procedure will be more or less the same, however, I will now utilize
 the function "round()" to account for MATLAB rounding issues as, my code's
 check will utilize a conditional "==" towards the end of this section. 
%}

% Calculate these first vector
vHD = randn(1,10)                                   ;

% Create a second vector and find the length of both vectors
uHD = randn(1,10)                                   ;
len_v_VecHD = sqrt(sum(vHD*vHD'))                   ;
len_u_VecHD = sqrt(sum(uHD*uHD'))                   ;

% Create a unit vector from vector u
unitVecHD = uHD/len_u_VecHD                         ;

% Find the projection line and round to six places
newLineMagnitudeHD = sum(vHD*unitVecHD')      ;
projectionUnitHD = newLineMagnitudeHD * unitVecHD   ;

% Find residual...
residualHD = vHD - projectionUnitHD                 ; 

%% Proving projection and unit lie on the same line.
%{
Because the dot product = the ||x|| * ||y|| * cos(theta_xy) where:
    x = a vector
    y = a vector
    theta_xy = the angle between them
We can find the angle between them by dividing the dot product by the
vector lengths and utilizing acosd()
%}

% Find the dot product of the projection onto the unit
dotUnitANDprojection = unitVecHD*projectionUnitHD'  ;

% Find these vectors lengths
length_unitVecHD = sqrt(unitVecHD*unitVecHD')       ;
length_projectionHD = sqrt(projectionUnitHD*projectionUnitHD')  ;

% Divide out the lengths and take the acos() of the result
cos_theta = dotUnitANDprojection/(length_unitVecHD*length_projectionHD) ;

% Round cos_theta to nearest 4 zeros to overcome rounding issues
cos_theta = round(cos_theta,4)  ;

% Display results based on the random vectors created...
if cos_theta == -1
    disp('The unit vector and the projection line are on the same line')
    disp('and they point in opposite directions.')
elseif cos_theta == 1
    disp('The unit vector and the projection line are on the same line')
    disp('and they point in opposite directions.')
end

%% FIX ME FIX ME FIX ME FIX ME FIX ME FIX ME FIX ME FIX ME FIX ME FIX ME FIX ME

% Show that the projection is orthogonal to the residual

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Find the dot product of the projection onto the unit
dotOrthogANDprojection = residualHD*projectionUnitHD' ; 

% Find these vectors lengths
length_orthogonal_VecHD = sqrt(residualHD*residualHD') ;
length_projectionHD = sqrt(projectionUnitHD*projectionUnitHD') ;

% Divide out the lengths and take the acos() of the result
cos_theta = dotOrthogANDprojection/(length_orthogonal_VecHD*length_projectionHD)

cos_theta < 10e-6

vecLenAngle(residualHD,projectionUnitHD) 


%% Make sure that the orthogonal vector plus the projection equals the original vector

%{ 
To show that the sum of the projection line and the orthogonal (residual)
line are equal to the original vector (vHD), we subtract vHD from these two
vectors sum. Then we utilize the function find() to determine where there
are zeros (because the conditional "==" will not work due to rounding
errors - we utilize the logical condition "< 10e-6").

This returns a list of the indices where the conditional requirement is
met.

As you can see, 

%}

difference_btw_vectors = (projectionUnitHD + residualHD) - vHD;
columnIndices = find(difference_btw_vectors < 10e-6) ;
equality_check = length(difference_btw_vectors) == length(columnIndices) ;

% Display results based on the random vectors created...
if equality_check == 1
    disp('The sum of the projection vector and the orthogonal vector equals the original vector')
else
    sprintf('The sum of the projection vector and the orthogonal vector DOES NOT equal the original vector\n it appears you have a problem with your calculations')
end


%% 3 Geometry of linear transformations

%{
a) Create the PlotVec2 function that takes in a matrix of height 2 and plots
each column vector from this matrix on a 2-dimensional axes. It sohuld
check that the matrix argument has height two, signaling an error if not.
Vectors should be plotted as a line from the origin to the vector position,
using circle or other symbol to denote the "head". It should also draw the
x and y axes, extending from -1 to 1. The two axes should be equal size, so
that horizontal units are equal to vertical units.

NOTE: Due to the question restrictions, with respect to the axes, it is
likely that some of the plotted vectors will reach past the border of the
plot figure space.
%}

% Lets give this function a matrix with the proper height...

rightSizeMAT = rand(2,12) ;
figure()
PlotVec2(rightSizeMAT)

% Lets transform the matrix and try again...
figure()
PlotVec2(incorrectSizeMAT')

% Now lets see what happens when you don't...
incorrectSizeMAT = rand(12,2) ;
PlotVec2(incorrectSizeMAT)

%% 3.b

%{
Write a second function vecLenAngle that takes two vectors as arguments and
returns the length of each vector, as well as the angle between them.
Decide how you would like to handle cases when one (or both) vectors have
zero length.

NOTE: 
    - This function can take row or column vectors.
    - 

%}

% First, we show that the function works with vectors of the same length...
vec1 = rand(12,1) ;
vec2 = rand(1,12) ;

vecLenAngle(vec1,vec2)

% Now, we show that the function throws errors for vectors of different length...

vec1 = rand(12,1) ;
vec2 = rand(1,10) ;

vecLenAngle(vec1,vec2)

% For good measure, we show that the function  works when the initial
% shapes match as well...

vec1 = rand(12,1) ;
vec2 = rand(12,1) ;

vecLenAngle(vec1,vec2)

%% 3.c 
%{
Generate a random 2x2 matrix M, and decompose it using the SVD, M = USV^T. Now
examine the action of this sequence of transformations on the two “standard basis”
vectors, {ˆe_1, ˆe_2}. Specifically, use vecLenAngle to examine the lengths and 
angle between two basis vectors ˆe_n, the two vectors V^T ˆe_n, the two vectors 
SV^T ˆe_n, and the two vectors USV^T ˆe_n. Do these values change, and if so, 
after which transformation? Verify this is consistent with their visual appearance 
by plotting each pair using plotVec2.
%}

randn

% Create a random 2x2 matrix
randMAT = randn(2) ;

subplot(2,3,1)
PlotVec2(randMAT);
title('Original Matrix');
grid on

% Decompes it utilizing the SVD
[U,S,V] = svd(randMAT) ;

% Create the two standard basis vectors
e_1 = [1;0] ;
e_2 = [0;1] ;

% Examine the length & angle between the basis vectors
subplot(2,3,2)
basis_values = vecLenAngle(e_1,e_2) ;
basis_combined = [e_1,e_2] ;
PlotVec2(basis_combined);
title('Basis Transform');
grid on

% Now do this at each step of the transformation...

% 1. V^T  : Examine the length & angle between the V^T*e_n

V_transform1 = V'*e_1 ;
V_transform2 = V'*e_2 ;
step1 = vecLenAngle(V_transform1,V_transform2)


subplot(2,3,4)
V_combined = [V_transform1,V_transform2] ;
PlotVec2(V_combined);
title('Step 1: V^T Transform');
grid on

% 2. SV^T : Examine the length & angle between the SV^T*e_n

SV_transform1 = S*V_transform1 ;
SV_transform2 = S*V_transform2 ;
step2 = vecLenAngle(SV_transform1,SV_transform2)

subplot(2,3,5)
SV_combined = [SV_transform1,SV_transform2] ;
PlotVec2(SV_combined);
title('Step 2: SV^T Transform');
grid on

% 3. USV^T: Examine the length & angle between the USV^T*e_n

USV_transform1 = U*SV_transform1 ;
USV_transform2 = U*SV_transform2 ;
step3 = vecLenAngle(USV_transform1,USV_transform2)

subplot(2,3,6)
USV_combined = [USV_transform1,USV_transform2] ;
PlotVec2(USV_combined);
title('Step 3: USV^T Transform');
grid on

test = round(USV_combined,6) == round(randMAT,6);

%{
As you can see from the lengths and angles, as well as the plots
themselves, the basis vectors get rotated onto the input 

%}

answer = 'Subjecting a matrix of basis vectors to the SVD\n of another vector returns the non-decomposed\n original vector.\n';

if test == true
    fprintf(answer)
else
    fprintf('Your code is broken, fix the calculations')
end

%% 3.d

%{
Generate a data matrix P with 65 columns containing 2-dimensional unit-vectors ˆu_n =
[cos(?n); sin(?n)], and ?n = 2pin/64, n = 0, 1, ... , 64.[Don’t use a for loop!
Create a vector containing the values of ?n. ] Plot a single blue curve through these
points, and a red star (asterisk) at the location of the first point. Consider the action
of the matrix M from the previous problem on this set of points. In particular, apply
the SVD transformations one at a time to full set of points (again, think of a way to do
this without using a for loop!), plot them, and describe what geometric changes you see
(and why).
%}

n = 1:64;
theta_func = @(n) (2*pi*n)/64;
thetas = theta_func(n);
unit_vec_func = @(theta) [cos(theta) ; sin(theta)];
COS_SINE = unit_vec_func(thetas);

plot(COS_SINE(1,:), COS_SINE(2,:), 'b-', 'Marker', '*', 'MarkerEdgeColor', 'r') ;
axis equal
title('Look ma'', I made a circle!')
xlabel('x-axis')
ylabel('y-axis')
legend('COS_ SINE Matrix Values')

%%

matrix = []

for i = 1:2
    value = (2*pi*i)/64 
    matrix = [matrix,[cos(value); sin(value)]] ;
   
end

x = [2,4,5]



plot(matrix(1,:), matrix(2,:), 'b-', 'Marker', '*', 'MarkerEdgeColor', 'r') ;
axis equal
title('Look ma'', I made a circle!')
xlabel('x-axis')
ylabel('y-axis')
legend('COS_ SINE Matrix Values')





%% 4 A Simple Visual Neuron

%{
Suppose a retinal neuron in a particular species of toad generates
responses that are a weighted sum of the (positive-valued) intensities of light
that is sensed at 6 localized regions of the retina. The weight vector is 
[1, 3, 8, 8, 3, 1]. 
%}

weighting_vector = [1;3;8;8;3;1] ;

%{ 

(a) Is this system linear? If so, how do you know? If not, provide a counter-example.

This system must be linear because the system's operations (a weighted sum) are
multiplication and addition. To be linear, a system must obey the rules
of superposition - homogeneity and additivity - both multiplication and
addition are allowed within a linear system.

We can also prove this by creating two different input vectors and testing
whether or not they obey the laws of super position.

%}

input_stimuli1 = [1,0,0,0,0,0] ;
input_stimuli2 = [0,1,0,0,0,0] ;


% Because input_stimuli2 is input_stimuli1 scaled by 2, putting these
% inputs into the system should return outputs which are also scaled by two
% - if the system is linear.


% Now we can check if output1 * 2 = output 2

if output1 * 2 == output2
    disp('The system is linear!')
else
    disp('The system is NOT linear!')
end


%{
(b) What unit-length stimulus vector (i.e., vector of light intensities) 
elicits the largest response in the neuron? 
- Explain how you arrived at your answer.
%}



%{
(c) What physically-realizable unit-length stimulus vector produces the smallest 
response in this neuron? 
- Explain. 
[hint: think about the geometry by visualizing a simpler version of the problem, in 2 dimensions]
%}

unit_weight_vec = weighting_vector / sqrt((weighting_vector'*weighting_vector));

unit_weight_vec' * weighting_vector

input_stimuli3 = [0,0,1,0,0,0] ;
input_stimuli4 = [0,0,0,1,0,0] ;
input_stimuli5 = [0,0,0,0,1,0] ;
input_stimuli6 = [0,0,0,0,0,1] ;

output1 = input_stimuli1 * weighting_vector
output2 = input_stimuli2 * weighting_vector
output3 = input_stimuli3 * weighting_vector
output4 = input_stimuli4 * weighting_vector
output5 = input_stimuli5 * weighting_vector
output6 = input_stimuli6 * weighting_vector


%% test

x = [2,1]
weights = [2,3]

frog = @(vec,weights) sum(vec*weights')
frog(x,weights)





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

**Extra credit: make your function recursive – instead of using a for
loop, have the function call itself, each time adding a new column to the matrix of previously
created orthogonal columns. To do this, you’ll probably need to write two functions (a
main function that initializes the problem, and a helper function that is called with a matrix
containing the current set of orthogonal columns and adds a new column until the number
of column equals the number of rows).

%}

% Create first vector
arb_vec = randn(1,10);

lengthVec = @(vec) sqrt(sum(vec*vec'))

% Get length
len_arb_vec = lengthVec(arb_vec)

% use that to normalize it
basis1 = arb_vec/len_arb_vec


%{
 You have to project onto EACH PREVIOUS BASIS, sum those up, and remove
 that total from the new random vector

%}

% Make second vector
arb_vec2 = randn(1,10);

% project onto the first basis
firstProjection = (arb_vec2*basis1')*basis1

% Find the orthogonal portion by removing the projection from the first
% vector
orthogonalVec = arb_vec2 - firstProjection

% Now take its length...
len_arb_vec2 = lengthVec(orthogonalVec)

% Normalize the new vector
basis2 = orthogonalVec/len_arb_vec2

% Check it is orthogonal to the last one
vecLenAngle(basis2,basis1)

% Make third vector
arb_vec3 = randn(1,10);

% project onto the first basis
secondProjection = (arb_vec3*basis2')*basis2
secondProjection2 = (arb_vec3*basis1')*basis1

% Find the orthogonal portion by removing the projection from the first
% vector
orthogonalVec = arb_vec3 - (secondProjection + secondProjection2)

% Now take its length...
len_arb_vec3 = lengthVec(orthogonalVec)

% Normalize the new vector
basis3 = orthogonalVec/len_arb_vec3

% Check it is orthogonal to the last one
vecLenAngle(basis1,basis2')

all_basis = [basis1 ; basis2 ; basis3]

figure 
hold on
for pp = 1:3
     plot3([0 all_basis(pp,1)],[0 all_basis(pp,2)],[0 all_basis(pp,3)])
end
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
view(3)
axis square
grid on
rotate3d on
legend('basis1','basis2','basis3')



%% Questions 6 - Null Space and Range Space

load('mtxExamples.mat')

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

find(mat2_zeros_mat < 10e-6) 
find(mat4_zeros_mat < 10e-6) 

%%

% Now we can generate a random vector based on the range space

% First we need to find the psuedo-inverse of S for each matrix...

mat_2_range = u2(:, mat2_range_index) ;
mat_4_range = u4(:, mat4_range_index) ;

S_2_diag = diag(s2)

S_2_inv_Vec = 1./S_2_diag(:)
S_2_inv = diag(S_2_inv_Vec)

x_input_vec = v2*(S_2_inv*(u2'*(randn*mat_2_range)))

mtx2*x_input_vec

S_4_diag = diag(s4)

S_4_inv_Vec = 1./S_4_diag(:)
S_4_inv = diag(S_4_inv_Vec)


random_mat2_vec_range = mat2_range_row * randn ;

result = 

mtx2* random_mat2_vec_range'

% Lets do this again for the fourth matrix...

mat_4_range
random_mat4_vec_range = mat4_range_row * randn ;

% Now we can transform them with U,S,V^transposed

mat2_range_mat = u2*(s2*(v2'*random_mat2_vec_range')) ;
mat4_range_mat = u4*(s4*(v4'*random_mat4_vec_range')) ;

u2

find(abs(mat2_range_mat) < 10e-6) ;
find(abs(mat4_range_mat) < 10e-6) ;


%% NOTES

%{
1) double check everything
2) change all of your "thing < 10e-6" to abs(thing < 10e-6"
3) make sure everything is complete.
4) double check all plots!
5) double check that the html file creates correctly!
6 Change the last question to be an addition of both orthogonal null basis
vectors
7 do the second part of the circle questions
%}


















