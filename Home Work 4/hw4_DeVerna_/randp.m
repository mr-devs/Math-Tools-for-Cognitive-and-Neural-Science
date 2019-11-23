function samples = randp(p,num)
%randp Generate samples from a given distribution

% Take the cumulative sum
cum = cumsum(p);

% Generate 'num' samples in the range (0,1), and then by the
% cumulative sum of the PDF. 
rand_vec = rand(1,num);

% Use histcounts to grab how many of each random number are in each bin specified 
[N,edges,samples] = histcounts(rand_vec,cum);

end


