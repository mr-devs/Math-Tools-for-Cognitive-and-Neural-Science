%% Math Tools I - HW #1 - Question 3: Geometry of Linear Transformations

% Purpose: This code has been written to wrestle with question 3 of the
% first homework of math tools I.

% Author: Matthew DeVerna

% Date: 09/24/19
%% (a) Create function PlotVec2

%{
a) Create the function PlotVec2 that takes in a matrix of height 2 and plots
each column vector from this matrix on a 2-dimensional axes. It sohuld
check that the matrix argument has height two, signaling an error if not.
Vectors should be plotted as a line from the origin to the vector position,
using circle or other symbol to denote the "head". It should also draw the
x and y axes, extending from -1 to 1. The two axes should be equal size, so
that horizontal units are equal to vertical units.

%}

% Due to the nature of the randn function, there can certainly be random
% numbers that are greater (or smaller) than 1 (or -1). In order to combat
% this aesthetic plotting issue. The absolute max value of each matrix will
% be found and then a square axis will be determined by that value + the
% below cushion value. 

cushion = .1 ;

% Lets give this function a matrix with the proper height...

rightSizeMAT = randn(2,12) ;
subplot(1,3,1)
PlotVec2(rightSizeMAT)
max_axis = max(max(abs(rightSizeMAT))) ;
cushion = .1 ;
xlim([-max_axis-cushion max_axis+cushion]) ;
ylim([-max_axis-cushion max_axis+cushion]) ;
title('First Try - Proper Input')

% Lets create a matrix of the wrong shape - i.e. height =/= 2...
incorrectSizeMAT = randn(12,2) ;

% We know this is wrong, but we can transform it and see that the plotting
% function still works...

subplot(1,3,2)
PlotVec2(incorrectSizeMAT')
max_axis = max(max(abs(incorrectSizeMAT))) ;
cushion = .1 ;
xlim([-max_axis-cushion max_axis+cushion]) ;
ylim([-max_axis-cushion max_axis+cushion]) ;
title('Second Try - Transposed Wrong Input')

% However, if we don't do that, we see how the function handles matrices
% that have been submitted to plot with the wrong dimensions.

% It seems like this error test keeps the rest of my code from running,
% feel free to uncomment and try it yourself, but all of my answers below
% do not publish if I keep this in here.

% subplot(1,3,3)
% title('Wrong Dimensions')
% PlotVec2(incorrectSizeMAT)

%% 3.b Create the function venLenAngle

%{
Write a second function vecLenAngle that takes two vectors as arguments and
returns the length of each vector, as well as the angle between them.
Decide how you would like to handle cases when one (or both) vectors have
zero length.

NOTE: 
    - This function can take row or column vectors.
    - Vectors with zero length are returned as such and it is assumed that
    the angle between the vectors is zero - as you cannot have an angle
    between two vectors when there is only one vector :).
%}

% First, we show that the function works with vectors of the same length...
vec1 = rand(12,1) ;
vec2 = rand(1,12) ;

[len1, len2, Angle] = vecLenAngle(vec1,vec2);

disp('Vectors with different dimensions')

fprintf('vecLen1 = %f \n vecLen2 = %f \n Angle = %f', len1, len2,Angle)

% For good measure, we show that the function works when the initial
% shapes match as well...

vec1 = rand(12,1) ;
vec2 = rand(12,1) ;

[len1, len2, Angle] = vecLenAngle(vec1,vec2);

fprintf('\n\nVectors with the same dimensions')

fprintf('vecLen1 = %f \n vecLen2 = %f \n Angle = %f', len1, len2,Angle)

% Now, we show that the function throws errors for vectors of different length...

% I commented this out as well, for the same reason as above. 

% vec1 = rand(12,1) ;
% vec2 = rand(1,10) ;
% 
% vecLenAngle(vec1,vec2)

%% 3.c Use PlotVec2 to put basis vectors through the svd transform step-by-step
%{
Generate a random 2x2 matrix M, and decompose it using the SVD, M = USV^T. Now
examine the action of this sequence of transformations on the two “standard basis”
vectors, {ˆe_1, ˆe_2}. Specifically, use vecLenAngle to examine the lengths and 
angle between two basis vectors ˆe_n, the two vectors V^T ˆe_n, the two vectors 
SV^T ˆe_n, and the two vectors USV^T ˆe_n. Do these values change, and if so, 
after which transformation? Verify this is consistent with their visual appearance 
by plotting each pair using plotVec2.
%}

% Create a random 2x2 matrix
randMAT = randn(2) ;

subplot(2,3,1)
PlotVec2(randMAT);
title('Original Matrix');
max_axis = max(max(abs(randMAT))) ;
cushion = .1 ;
xlim([-max_axis-cushion max_axis+cushion]) ;
ylim([-max_axis-cushion max_axis+cushion]) ;
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
max_axis = max(max(abs(basis_combined))) ;
cushion = .1 ;
xlim([-max_axis-cushion max_axis+cushion]) ;
ylim([-max_axis-cushion max_axis+cushion]) ;
title('Basis Transform');
grid on

% Now do this at each step of the transformation...

% 1. V^T  : Examine the length & angle between the V^T*e_n

V_transform1 = V'*e_1 ;
V_transform2 = V'*e_2 ;

[len1, len2, Angle] = vecLenAngle(V_transform1,V_transform2);

fprintf('\n\nStep 1: V'' Transform  (Rotate)\n')

fprintf('vecLen1 = %f \n vecLen2 = %f \n Angle = %f', len1, len2,Angle)

subplot(2,3,4)
V_combined = [V_transform1,V_transform2] ;
PlotVec2(V_combined);
max_axis = max(max(abs(V_combined))) ;
cushion = .1 ;
xlim([-max_axis-cushion max_axis+cushion]) ;
ylim([-max_axis-cushion max_axis+cushion]) ;
title('Step 1: V'' Transform');
grid on

% 2. SV^T : Examine the length & angle between the SV^T*e_n

SV_transform1 = S*V_transform1 ;
SV_transform2 = S*V_transform2 ;
[len1, len2, Angle] = vecLenAngle(SV_transform1,SV_transform2);

fprintf('\n\nStep 2: SV'' Transform (Stretch/Shrink)\n')

fprintf('vecLen1 = %f \n vecLen2 = %f \n Angle = %f', len1, len2,Angle)

subplot(2,3,5)
SV_combined = [SV_transform1,SV_transform2] ;
PlotVec2(SV_combined);
max_axis = max(max(abs(SV_combined))) ;
cushion = .1 ;
xlim([-max_axis-cushion max_axis+cushion]) ;
ylim([-max_axis-cushion max_axis+cushion]) ;
title('Step 2: SV'' Transform');
grid on

% 3. USV^T: Examine the length & angle between the USV^T*e_n

USV_transform1 = U*SV_transform1 ;
USV_transform2 = U*SV_transform2 ;
[len1, len2, Angle] = vecLenAngle(USV_transform1,USV_transform2);

fprintf('\n\nStep 3: USV'' Transform (Rotate)\n')

fprintf('vecLen1 = %f \n vecLen2 = %f \n Angle = %f', len1, len2,Angle)


subplot(2,3,6)
USV_combined = [USV_transform1,USV_transform2] ;
PlotVec2(USV_combined);
max_axis = max(max(abs(USV_combined))) ;
cushion = .1 ;
xlim([-max_axis-cushion max_axis+cushion]) ;
ylim([-max_axis-cushion max_axis+cushion]) ;
title('Step 3: USV'' Transform');
grid on

test = round(USV_combined,6) == round(randMAT,6);

answer = '\n\nSubjecting a matrix of basis vectors to the SVD\n of another vector returns the non-decomposed\n original matrix.\n';

if test == true
    fprintf(answer)
else
    fprintf('Your code is broken, fix the calculations')
end

%% 3.d Create a Circle Matrix, Plot it, Then Plot Again and Again after SVD Transforms from Previous Problem Matrix

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

n = 0:64;
thetas_ = (2*pi) * n ;
thetas = thetas_ / 64 ;
COS_SINE = [cos(thetas) ; sin(thetas)] ;

figure ;
subplot(2,3,2) ;
plot(COS_SINE(1,1), COS_SINE(2,1), 'Marker', '*', 'MarkerEdgeColor', 'r', 'MarkerSize', 10)
hold on
plot(COS_SINE(1,:), COS_SINE(2,:), 'b-') ;
axis equal
title('Orig. COS_ SINE Matrix Values')
max_axis = max(max(abs(COS_SINE))) ;
cushion = .1 ;
xlim([-max_axis-cushion max_axis+cushion]) ;
ylim([-max_axis-cushion max_axis+cushion]) ;
xlabel('x-axis')
ylabel('y-axis')

% Plot the matrix step by step after each SVD transformation step

V_transform_circle = V'*COS_SINE ;

subplot(2,3,4) ;
plot(V_transform_circle(1,1), V_transform_circle(2,1), 'Marker', '*', 'MarkerEdgeColor', 'r', 'MarkerSize', 10)
hold on
plot(V_transform_circle(1,:), V_transform_circle(2,:), 'b-') ;
axis equal
title('V'' Transform')
max_axis = max(max(abs(V_transform_circle))) ;
cushion = .1 ;
xlim([-max_axis-cushion max_axis+cushion]) ;
ylim([-max_axis-cushion max_axis+cushion]) ;
xlabel('x-axis')
ylabel('y-axis')

SV_transform_circle = S*(V'*COS_SINE) ;

subplot(2,3,5) ;
plot(SV_transform_circle(1,1), SV_transform_circle(2,1), 'Marker', '*', 'MarkerEdgeColor', 'r', 'MarkerSize', 10)
hold on
plot(SV_transform_circle(1,:), SV_transform_circle(2,:), 'b-') ;
axis equal
title('SV'' Transform')
xlabel('x-axis')
ylabel('y-axis')
max_axis = max(max(abs(SV_transform_circle))) ;
cushion = .1 ;
xlim([-max_axis-cushion max_axis+cushion]) ;
ylim([-max_axis-cushion max_axis+cushion]) ;

USV_transform_circle = U*(S*(V'*COS_SINE)) ;

subplot(2,3,6) ;
plot(USV_transform_circle(1,1), USV_transform_circle(2,1), 'Marker', '*', 'MarkerEdgeColor', 'r', 'MarkerSize', 10)
hold on
plot(USV_transform_circle(1,:), USV_transform_circle(2,:), 'b-') ;
axis equal
title('USV'' Transform')
xlabel('x-axis')
ylabel('y-axis')
max_axis = max(max(abs(USV_transform_circle))) ;
cushion = .1 ;
xlim([-max_axis-cushion max_axis+cushion]) ;
ylim([-max_axis-cushion max_axis+cushion]) ;

answer = 'First Step = Rotate \n\nSecond Step = Stretch \n\nThird Step = Rotate Again';

sprintf(answer)

