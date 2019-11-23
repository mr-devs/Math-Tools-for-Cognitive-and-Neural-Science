%% Math Tools I - HW #1 - Question 1

% Purpose: This code has been written to wrestle with question 1 of the
% first homework of math tools I.

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

answer = [('As we can see, this system does not follow the laws of superposition'),...
    (' (specifically homogeneity). Thus, we conclude that this system is not linear.')]

%% Question 1: System 2:

%{
Since the scope of this question's programming complexity is quite simple,
I will use the same variables as we move through the rest of this problem. 
The linear nature of this script allows me to do so without worrying about 
workspace/variable issues - and offers the benefit of simple, clear and
understandable variables.
%}

input1 = [6;3]      ;
input2 = [-2;-1]    ;
output1 = [12;12]   ;
output2 = [-6;-6]   ;

% Because scaling input1 by -1/3 gives input2 exactly, the same affect
% should be seen for the outputs - if the system is linear.

if input1 * (-1/3) == input2                        % Set conditional if-statement
    disp('Input 1 * -1/3 equals Input 2')           % If the condition is accurate, print this.

else
    disp('Input 1 * -1/3 DOES NOT equal Input 2')   % If the condition is NOT accurate, print this.
end

% Lets try the same for the outputs...

if output1 * (-1/3) == output2         % Set conditional if-statement
    disp('Output 1 * -1/3 equals Output 2')          % If the condition is accurate, print this.
else
    disp('Output 1 * -1/3 DOES NOT equal Output 2')  % If the condition is NOT accurate, print this.
end

answer = [('Once again, we can see that this system does not follow the laws of superposition'),...
    (' (specifically homogeneity). Thus, we conclude that this system is not linear')]

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
    disp('Input3 * 2/3 - Input1 equals Input 2')            % If the condition is accurate, print this.

else
    disp('Input3 * 2/3 - Input1 DOES NOT equals Input 2')   % If the condition is NOT accurate, print this.
end

% Lets try the same for the outputs...

if (output3 * 2/3) - output1 == output2                         % Set conditional if-statement
    disp('Output 3 * 2/3 - Output1 equals Output 2')            % If the condition is accurate, print this.
else
    disp('Output 3 * 2/3 - Output1 DOES NOT equals Output 2')   % If the condition is NOT accurate, print this.
end

answer = [('Daaamn, my (wo)man, these systems just DO NOT want to follow the laws of superposition'),...
    (' (specifically homogeneity). Thus, we conclude that this system is not linear.')]

%% Question 1: System 4

input1 = [2 ; 4];
input2 = [-2; 1];
output1 = 0 ;
output2 = 3 ;

%{ 
We can determine the values of the 2-by-1 vector that is the linear system.
The shape of the system is determined by the shape of the input and output.

Input   -->  System -->  Output
(1 x 2) --> (2 x 1) -- (1 x 1)
--> Take the "outer elements" of the matrix shapes to form the end result.

We can determine this vector's values through the elimination method of
solving linear equations:

[2 ; 4] -->  0 
[-2; 1] -->  3    can be written as:

2x + 4y = 0    and
-2x + y = 3
|--> Because the linear system's effect is equivalent to the dot product.

Then, combining these equations into one gives you:

5y = 3

Solving for y gives you --> y = .6 or (3/5)

Plug in .6 to one of the equations and solve for x...

2x + 4(3/5) = 0
2x + 2.4    = 0
2x          = -2.4
x           = -1.2

Based on this solution, it looks as though these input/output pairs create
a linear system which should successfully operate on our pairs in a
predictable manner - such that we get the predicted output pairs.

Lets see if we can test this below in the same manner as earlier problems...
%}

% Create our test vector...
testVec = [-1.2,.6] ;


if (testVec * input1) == output1                                            % Set conditional if-statement
    disp('The dot product of testVec and Input1 equals Output1')            % If the condition is accurate, print this.

else
    disp('The dot product of testVec and Input1 DOES NOT equals Output1')   % If the condition is NOT accurate, print this.
end

% Lets try the same for the input2 and output2...

if (testVec * input2) == output2                                            % Set conditional if-statement
    disp('The dot product of testVec and Input2 equals Output2')            % If the condition is accurate, print this.

else
    disp('The dot product of testVec and Input2 DOES NOT equals Output2')   % If the condition is NOT accurate, print this.
end

answer = [('By testing the matrix that we found, we can see this system is linear')]

% Furthermore, we can scale these inputs and feed them into the system to
% further verify them.

if (testVec * (3*input1)) == 3*output1                                            % Set conditional if-statement
    disp('The dot product of testVec and Input1 (scaled by 3) equals Output1')            % If the condition is accurate, print this.

else
    disp('The dot product of testVec and Input1 (scaled by 3) DOES NOT equals Output1')   % If the condition is NOT accurate, print this.
end

% Lets try the same for the input2 and output2...

if (testVec * (3*input2)) == 3*output2                                            % Set conditional if-statement
    disp('The dot product of testVec and Input2 (scaled by 3) equals Output2')            % If the condition is accurate, print this.

else
    disp('The dot product of testVec and Input2 (scaled by 3) DOES NOT equals Output2')   % If the condition is NOT accurate, print this.
end

answer_check = [('Even after scaling the inputs and outputs, and putting them through'),...
    (' the linear system that we found - this system still holds up to the laws of superposition.')]


%% Question 1: System 5

input1 = 0;
ouput1 = [1;2];

answer = sprintf([('We can tell right away that this system is not linear because we are starting at'),...
    (' the origin. By definition,\n there is no linear combination that can take us away'), ...
    (' from the origin and/or add a new dimension')])


