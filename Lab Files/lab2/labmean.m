function avg = labmean(inputVec)
% A custom implementation of mean(). Takes in a vector of values and
% returns the average.

% Define an accumulator variable to serve as a running count of the values in the vector
currSum = 0;
% Loop to take the sum of the vector values
for ii = 1:length(inputVec)
    currSum = currSum + inputVec(ii);
end


avg = currSum/length(inputVec);

end