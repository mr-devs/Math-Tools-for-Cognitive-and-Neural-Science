function sum = computeSum(start, num, space)
% Takes in 3 inputs: start, num, and space. Beginning at start, this 
% returns the sum of num numbers, each at an interval of space from the 
% previous one.

% Set two variables: one to hold the current number in the series and one 
% to keep track of the cumulative sum
curr_num = start;
sum = start;

% Loop through num numbers
for i = 1:num-1
    curr_num = curr_num + space;
    sum = sum + curr_num;
end
