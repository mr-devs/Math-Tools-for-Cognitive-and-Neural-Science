%% MATH TOOLS 2019 LAB 2: MATLAB Intro II

%% 1 IF/ELSE statements

% Used to execute a code chunk only when certain conditions are met

% if condition
%   some code
% elseif condition
%   more code
% else
%   still more code
% end

%%
% * Example: day of the week

% First let's get the day of the week
[day_num, day_name] = weekday(datetime)

% Now we'll check if we have lab or lecture today
if strcmp(day_name, 'Fri') == 1    % if today is Friday, then we have lab
    class_type = 'lab';
else                               % on any other day, we have lecture
    class_type = 'lecture';
end

% Display the results as a string
disp(['Today''s class is: ' class_type])

% Note: this is the first time we've worked with strings in this class.
% These differ from integers or vectors in many way. For example, we define
% strings with '' or "" and comparison between strings is done with the
% function strcmp() rather than "==".

% We can modify our code further to be more specific, since we actually
% only have lecture on Tuesdays and Thursdays
if strcmp(day_name, 'Fri') == 1                         % if today is Friday, then we have lab
    class_type = 'lab';
elseif strcmp(day_name, 'Tue') || strcmp(day_name, 'Thu') % if today is Tuesday OR Thursday
    class_type = 'lecture';
else                               
    class_type = 'none';                                % on any other day, we don't have class
end

disp(['Today''s class is: ' class_type])


% Another note: use "||" to mean OR and "&&" to mean AND in conditional
% statements.

% Now test our code fully by manually setting day_name ('Mon', 'Tue', 'Wed', 'Thu', 'Fri')
day_name = 'Wed';
if strcmp(day_name, 'Fri') == 1                         % if today is Friday, then we have lab
    class_type = 'lab';
elseif strcmp(dayName, 'Tue') || strcmp(dayName, 'Thu') % if today is Tuesday OR Thursday
    class_type = 'lecture';
else                               
    class_type = 'none';                                % on any other day, we don't have class
end

disp(['Today''s class is: ' class_type])

%%
% * Exercise: boundary conditions

% Write code that generates a random integer between -50 and 50 (hint: use
% randi()) and displays it. Then, convert the given number into the range
% from -10 to 10 and display it (i.e. 25 becomes 10, -15 becomes -10,
% etc.)

%% 2 WHILE loops

% Used to run a section of code to be run for a certain amount of
% time, or repeatedly as long as a certain condition is met.

% while condition
%   some code
% end

%%
% * Example: timing code

% To see how long something takes, we can measure elapsed time with "tic"
% and "toc". These are useful for estimating code efficiency as well.

tic             % start time
while toc < 2   % while the elapsed time is less than 2 seconds
    disp(toc)   % display elapsed time
end

% Note: you can clear the command window  with "clc"

%%
% * Example: randomly distributed numbers

% While loops are also helpful when you're not sure a priori how many
% iterations you need

% Draw 2 uniformly distributed random numbers in the interval (0,1)
sample1 = rand
sample2 = rand

% Let's say we want sample2 to be smaller than sample1. We can use a while
% loop to repeatedly pull samples until the condition is met.

while sample2 >= sample1
    sample1 = rand
    sample2 = rand
end

% Now let's incrementally add to sample2 until its larger than sample1
while sample2 < sample1
    sample2 = sample2 + .01
end

% An alternate way to do this is to set while to 1 and have an if statement
% in the while loop. Let's try  that.

% First reset the numbers:
sample1 = rand; 
sample2 = sample1/2; % ensure sample2 starts as smaller than sample1

while 1
    if sample2 > sample1
        disp([num2str(sample2) ' is larger than ' num2str(sample1)])
        % num2str() converts numbers into strings so that they can be used
        % by functions that require strings
        break
    end
    sample2 = sample2 + .01
end

% Finally, let's say we did this problem but accidentally subtracted from
% sample2 rather than added. What happens?

sample1 = rand; 
sample2 = sample1/2;

while sample2 < sample1
    sample2 = sample2 - .01
end

% This creates an infinite loop. Type "ctrl+c" in the command window
% to break.

%%
% * Exercise: cumulative random integers

% Write some code to generate random integers between 1 and 20 until they
% add up to more than 50. At every iteration through
% of your while loop, display the current sum of the integers.


%% 3 FOR loops

% Used when you know how many times you want to iterate. By default, each 
% iteration through the for loop updates the counter by 1.

% for variable = start:stride:stop
%   some code
% end

%% 
% * Example: counting

% Increasing from 1 to 5
for ii = 1:5
    disp(ii)
end

% Decreasing from 5 to -2
for ii = 5:-1:-2 % decrease
    disp(ii)
end

% Increasing by non-integers
for ii = 0:0.1:1
    disp(ii)
end

%% 
% * Example: fibonacci sequence
% This is a function where the next value in a sequence depends on the
% previous values.

% Initialize your matrix. Why do we do this?
fib = [];
% Iterate 10 times
for ii = 1:10
    if length(fib) == 0       % if the current length of fib is 0
        fib(ii) = 0;
    elseif length(fib) == 1   % if the current length of fib is 1
        fib(ii) =  1;
    else
        fib(ii) = fib(ii-1) + fib(ii-2);
    end
end

% Other ways to initialize a matrix in matlab: zeros() and ones()

% Now let's repeat this but step through the loop incrementally using 
% breakpoints to see what's going on. This is a useful way to debug more
% complicated code.
fib = [];
for ii = 1:10
    if length(fib) == 0       % if the current length of fib is 0
        fib(ii) = 0;
    elseif length(fib) == 1   % if the current length of fib is 1
        fib(ii) =  1;
    else
        fib(ii) = fib(ii-1) + fib(ii-2);
    end
end

%% 
% * Example: simple data vectors

% Let's create two data vectors, each with 10 elements and store those as a
% single matrix
nVec = 2;
vecLength = 10;
vecs = rand(vecLength, nVec);

% Display the vectors as a stem plot
figure; hold on;
% subplot 1
subplot(2,1,1)
stem(vecs(:,1))
title('subplot 1')
% subplot 2
subplot(2,1,2)
stem(vecs(:,2))
title('subplot 2')

%% 
% * Exercise: more data vectors

% Write code to display 20 vectors instead of just 2. Have each of the
% subplots be titled like above, and have more than 1 column of subplots.


%% 4 Indexing

% Loops are useful and often necessary, but they aren't very efficient in 
% matlab. When possible, use indexing to change multiple values in a vector 
% at the same time.

% Reminder: you can index into a vector by using parentheses
x = rand(5,1)   % create a 5-vector
x(3)            % select the 3rd index
x(3) = x(4)     % set the value at the 3rd index to the value at the 4th index
x(3) == x(4)    % evaluates to 1

%%
% * Exercise: indexing vs for loops

% Using a for loop, create a vector that is 1e6 (1e6 = 1*10^6 = 1000000)
% rows where all the values are zero except for all indices that are
% multiples of 3, where the value should be 1. Then do the same thing 
% using indexing, and use tic/toc o show how much faster indexing is.

%% 5 Functions

% We've been using functions throughout these labs: dot(), sum(), rand(),
% etc. are all native to matlab. If there is a computation that you
% will need to use repeatedly, it is good practice to define your own
% function.

% function output = function_name(input)
%   some code
% end

% Some basics:
% - to create a new function, go to the EDITOR tab above and click 
%   new -> function
% - the inside of the function will flexibly respond to whatever the input
%   argument is
% - the intermediate values will not exit the function
% - only what you set equal to the output will be output
% - save the file with the same name as the function definition

%%
% * Example: mean function

% To keep things basic, lets re-create the mean() function. Note that you 
% need to name it something other than mean() so you don't overwrite
% matlab's native implementation. We'll verify that this works by comparing
% output from our function and matlab's.

rand_vec = rand(1,5);
avg = labmean(rand_vec)
avg == mean(rand_vec)

%%
% * Exercise: 

% Write a function, computeSum, that takes in 3 inputs: start, num, and
% space. Beginning at start, this should return the sum of num numbers
% that are space apart. Test your code on the following computations.

% The sum of the numbers 1 through 10 is 55
computeSum(1,10,1)

% The sum of the first 5 even numbers is 30
computeSum(2,5,2)

% The sum of the number -1 through -10 is -55
computeSum(0,11,-1)

%% 6 Review of vectors and matrices

% Inner products
rowVector = [1 2 3 4];
columnVector = [5;6;7;8];               % or [5 6 7 8]'
innerProduct = rowVector*columnVector
dot(rowVector,columnVector)             % same
sum(rowVector.*columnVector')           % same

% Matrix multiplication 
A = [1 2; 3 4];
B = [5 6; 7 8];
C = A*B
D = B*A

% Matrix multiplication does not commute. Neither do dot products, as a
% special case of matrix multiplication. Only element-wise multiplication
% commutes.

E1 = A.*B
E2 = B.*A
E1 == E2

%%
% * Exercise: combining a while loop, if statements, and dot products

% Let's see how many iterations it takes us to find 10 orthogonal vectors 
% by random trial and error, and keep track of what those pairs are. 
% Replace the provided pseudocode with functioning matlab, which should
% then work with the plotting code at the bottom of the section to
% visualize those vectors.

% Set counters to keep track of the number of loop iterations and number
% of orthogonal vectors. Also initialize two matrices to keep track of your
% row and column vectors. 

row_container = zeros(10,2);
column_container = zeros(2,10);

num_ortho = 0;
num_iter = 0;

while num_ortho < 10
    % Draw a 2D row vector randomly from a normal distribution (hint: use
    % randn())
    % Draw a 2D column vector randomly from a normal distribution
    % Increment loop counter
    
    % Test if the random vectors are orthogonal (hint: test for an inner
    % product very close to 0 rather than 0, because matlab has some
    % numeric variability). If so, increase the counter for the number of 
    % orthogonal vectors, and store each of the two vectors in the
    % corresponding index of their array.
end

% Let's test if the above code worked. First, print the number of
% iterations it took to find the 10 pairs.
disp(['It took ', num2str(num_iter), ' iterations to find 10 orthogonal vectors']);

% Code for plotting the orthogonal pairs
figure; hold on;
for ii = 1:10
    subplot(2,5,ii)
    plot([0 row_container(ii,1)],[0 row_container(ii,2)],'color','b')
    axis([-2 2 -2 2])  
    
    hold on; % put hold on, otherwise the second vector will overwrite the first one
    plot([0 column_container(1,ii)],[0 column_container(2,ii)],'color','r')
    hold off;
    title(['pair ' num2str(ii)])
end
