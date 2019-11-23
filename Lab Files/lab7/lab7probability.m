%% MATH TOOLS LAB 7: PROBABILITY THEORY

%% 1 Summary statistics

% Let's say we have a single data point x_i = 1. We're going to
% plot the l0, l1, l2 norm discrepancies to a statistic s.

% Define a vector s from -1 to 3 with 41 points (each 0.1 apart)
x_i = 1;
s = linspace(-1,3,41);

% L0 norm: index where x_1 occurs in s
ind = find(s==x_i);
l0norm = ones(1,41);
l0norm(ind) = 0;
plot(s, l0norm, 'LineWidth', 2)
hold on;

% L1 norm: absolute value of the difference between x_i and s
l1norm = abs(x_i - s);
plot(s, l1norm, 'LineWidth', 2)

% L2 norm: squared difference between x_i and s
l2norm = (x_i - s).^2;
plot(s,l2norm, 'LineWidth', 2)

% Plot all together
legend('l0', 'l1', 'l2')
xlabel('Summary statistic value')
ylabel('Discrepancy to data point x_i = 1')

% In this case the statistic s = x_1. This is obvious but the shape of the
% curves gives you some intuition about what each norm statistic means. Now
% we can think about this in aggregate over a bunch of data points: this 
% leads us to mean, median, and mode.

%% Exercise: norms

% Say x = [0.5,1,1,1,3]. Compute the aggregate norms and plot for
% s = linspace(-1,3,41). Then, use min to find the minimum values of each
% norm, and compare those to the summary statistic they represent. Think
% about which summary statistics are good representations of the data and
% why that is.



%% 1 Summary statistics, cont.
% We can see this again in a real example comparing a Gaussian to
% a skewed distribution (using randn and pearsrnd)

figure; hold on;
x = randn(10000,1);
histogram(x)
meanval = mean(x);
medianval = median(x);
line([meanval, meanval], ylim, 'LineWidth', 2, 'Color', 'r');
line([medianval, medianval], ylim, 'LineWidth', 2, 'Color', 'c');
legend('datahist','mean','median')

figure; hold on;
x = pearsrnd(0,1,0.8,3,10000,1);
histogram(x)
meanval = mean(x);
medianval = median(x);
line([meanval, meanval], ylim, 'LineWidth', 2, 'Color', 'r');
line([medianval, medianval], ylim, 'LineWidth', 2, 'Color', 'c');
legend('datahist','mean','median')

% The median here is a better statistic for the skewed distribution even
% though for a normal distribution the mean and median are both equally
% good.


%% 2 Distributions and sampling

% Take a die that has a uniform probability of showing one of the 6 sides
% What is the distribution of X?  Uniform

figure; hold on;
x = 1:1:6;
y = [1/6,1/6,1/6,1/6,1/6,1/6];
plot(x,y,'LineWidth',2);
title('Distribution for a dice roll')
xlabel('x')
ylabel('P(X=x)')

% Given the pdf of a function, how do we generate samples?

% Option #1: randsample(n,k) generates k samples from 1:n
% Generate 50 samples and histogram the output 
figure; hold on;
samps = randsample(6,50,1);
histogram(samps);
title('Sampling for a dice roll (randsample)')
xlabel('Value')
ylabel('Frequency')

% Option #2: generalize to any pdf
% Get the cumulative sum
c = cumsum(y);

% Generate 50 numbers in the range (0,1), and then use histcounts to grab
% how many of each random number are in each bin specified by the
% cumulative sum of the PDF.
rand_vec = rand(1,50);
[N,edges,samples] = histcounts(rand_vec,c);

figure; hold on;
histogram(samples);
title('Sampling for a dice roll (generalized)')
xlabel('Value')
ylabel('Frequency')

%% Exercise: Mean, expected value, and variance

% Calculate the expected value of rolling a single dice. Then, calculate
% the empirical mean from the samples you generated above. Are these equal?
% Do the same for the variance and the sample variance.  Then, compute the 
% average of the sample means and the average of the sample variances over 
% 100,000 different generations of samples.



%% 3 Covariance matrix

% Generate random 2D points with mean at (0,0) and variance 1, also called
% white noise.
figure; hold on;
data = normrnd(0,1,[500,2]);
scatter(data(:,1),data(:,2))
title('Generated data');
xlabel('x');
ylabel('y');

% This has the identity matrix as its covariance matrix, but let's check
% that that's true.

x_mean = mean(data(:,1));
y_mean = mean(data(:,2));
cov_mat11 = mean((data(:,1)-x_mean).*(data(:,1)-x_mean));
cov_mat12 = mean((data(:,1)-x_mean).*(data(:,2)-y_mean));
cov_mat21 = mean((data(:,2)-y_mean).*(data(:,1)-x_mean));
cov_mat22 = mean((data(:,2)-y_mean).*(data(:,2)-y_mean));
cov_mat = [cov_mat11 cov_mat12 ; cov_mat21 cov_mat22]


%% Exercise: Linear transformations

% We want to characterize how transformations affect both the data and the
% covariance matrix. First, use the scaling matrix provided to stretch the
% x and y components of the data, plot it, and then re-calculate the covariance
% matrix (you can just use cov()). How did the matrix change? You should 
% be able to reason about this.

% Then, repeat this after applying the provided rotation matrix to your scaled data.
% Now it will be hard to decompose the covariance matrix back into rotation
% scaling matrices, and in fact you can't do this "by eye". Instead,
% decorrelate the data using its eigendecomposition with the code provided
% and try to understand why this works.

scaling_matrix = [0.7 0; 0 3.4];
theta = 0.77*pi;
rotation_matrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];




%% PCA whitening
% Removing correlations and normalizing features is an important piece of
% many learning algorithms and is something that the visual system must do
% because we would like to remove redundant information and make each
% feature as informative as possible. This also makes data such that the
% features all have unit variance.

% How does whitening the data work? 
% We have our data X and its covariance matrix, Cx. Our goal is to convert
% it into the dataset Y, which has identity covariance. In other words, we
% want to change the spherical, rotated shape of the data back into a
% spherical, oriented point cloud. The derivation to do this is as follows:

% Start with our equation:
% Y = MX

% Replace Y with the identity matrix, and X with its covariance matrix
% multiplied by the transpose of M (remember, we want to turn our cov
% matrix into the identity via a transform with M):
% I = MCxMt

% Now we eigen-decompose Cx into VDVt, where D is a diagonal matrix of
% eigenvalues (variances) and V are the eigenvevtors:
% I = MVDVtM

% D is our variance and D^1/2 is the standard deviation of the matrix, so
% we can replace D:
% I = MV(D^1/2)(D^1/2)VtMt

% Now we have our solution, as plugging this M into the equation causes
% everything to cancel and gives us the identity matrix. So, when the
% matrix M is multiplied by X, our data is "whitened" (turns from an
% ellipsoid to a sphere).
% M = (D^(-1/2)*V^T)

% Here's the implementation:
% Re-rotate the data to decorrelate it using the PCA axes. Make sure the
% variables rcov_mat and rotated_data match your final covariance matrix
% and rotated data from the exercise above.
[V,D] = eig(rcov_mat);
oriented_data = V*rotated_data';
figure; hold on;
plot(oriented_data(1,:), oriented_data(2,:), '+')
xlim([-4,4])
ylim([-10,10])
title('Re-oriented data');
xlabel('x');
ylabel('y');

ocov_mat = cov(oriented_data')

% Scale the variances of each feature to be the same
Dsqinv = [1/sqrt(D(1,1)),0 ; 0, 1/sqrt(D(2,2))];
whitened_data = Dsqinv*oriented_data;
figure; hold on;
plot(whitened_data(1,:), whitened_data(2,:), '+')
title('Uncorrelated data');
xlabel('x');
ylabel('y');

wcov_mat = cov(whitened_data')

% Another note: In one of the homework #4 problems, we start from an
% identity covariance matrix and turn it into an ellipsoid-shaped
% distribution with covariance Cx. Here, we're going in the reverse
% direction.