clear
close all

%% 1-generate 2-D gaussian distribution
% given a correlation coefficient, how to generate data samples? 
target_r = .75;
nsample = 1000;
%% way 1: generate distribution from covariance matrix using matlab function
% -- for homework of course you shouldn't use that. 
std1 = 1;
std2 = 1;
std12 = std1*std2*target_r;
covmat = [std1,std12;std12,std2];% "semi-definite" matrix

meanvec = [0,0];
R = mvnrnd(meanvec,covmat,nsample);
figure;
plot(R(:,1),R(:,2),'.')
title(['r^2 =' num2str(corr2(R(:,1),R(:,2)))])

%% way 2: generate data of unit variance then transform
% first, generate two orthogonal dataset that their covariance is identity
X1 = randn(nsample,1);
X2 = randn(nsample,1);
R0 = [X1,X2]; % or, just use randn(nsample,2)
% check: is X1 and X2 independent?
disp(R0'*R0/(nsample-1)) %way1: check covariance close to I
disp(corr2(X1,X2))%way2: check correlation coefficient close to 0
[u,s,v]=svd(R0);
disp(s(1:2,:)/sqrt(nsample)) %way3: check s matrix of original data svd close to I
if abs(corr2(X1,X2))<.01
    disp('yes the covariance matrix of R0 is identity')
end

% want the new variance to be sth new: C = M'IM. What should be M?
% C = USU' 
% C = Y'Y = (XM)'XM = M'IM = M'M
% => M = sqrt(S)U'.

[u,d,v] = svd(covmat);
M = d.^.5 * v';
R = R0*M;

figure;
plot(R(:,1),R(:,2),'.')
title(['r^2 =' num2str(corr2(R(:,1),R(:,2)))])

%% way 3: use the geometric intuition of coefficience
X1 = randn(nsample,1);
X2 = randn(nsample,1);
Y1 = X1+X2*tan(acos(target_r));
R = [X1,Y1];
figure;
plot(R(:,1),R(:,2),'.')
title(['r^2 =' num2str(corr2(R(:,1),R(:,2)))])

disp([norm(X1),norm(X2)])
disp(X1'*X2/norm(X1)/norm(X2)) 

% FYI: can you use this method to generate data given the full coefficient matrix?


%% 2-operation of random variables
% What is the mean and variance of:
% Y = X1+X2
% where X1,X2 are independently drawn from 1-D gaussian
%% what's your intuition? should it become bigger or smaller?
mean1 = 0;
std1 = 5;
mean2 = 20;
std2 = 2;
X1 = std1*randn(nsample,1)+mean1; % generate normal variables of given parameters
X2 = std2*randn(nsample,1)+mean2;
%% way1: linear transformation
% E(Y) = E(X)M
% E(Y^2) = ME(X^2)M'
M = [1;1];
meanY = [mean1,mean2]*M;
varY = M'*[std1^2,0;0,std2^2]*M;
sprintf('mean of X1+X2:%.2f, variance: %.2f',meanY,varY)

%% way2: vector perspective for determining variance
sprintf('variance: %.2f',std1^2 + std2^2)
%% way3: do simulation
Y1 = X1+X2;
sprintf('mean of X1+X2:%.2f, variance: %.2f,std: %.2f',mean(Y1),var(Y1),std(Y1))

%% 3-conditional, marginal probability (prosecutor fallacy)
% see slides.
p_table = 196/10194
p_bayes = .98*.0002/.0102

%% 4-Bayesian belief update
%% 4.1-binomial
N = 12; % data from the class
p = 5/12;
Bino = makedist('Binomial','n',N,'p',p);
sprintf('mean=%.2f, standard dv = %.2f',mean(Bino),std(Bino))

%% 4.2-Beta and Bernoulli distribution
expBern = p;
varBern = p*(1-p);
sprintf('mean=%.2f, standard dv = %.2f',N*expBern,sqrt(N*varBern))

%% 4.3-Inferring posterior of the bias rate
Data1 = [0,0,0,0,0,0,0,1,1,1,1,1]; % data from the first sample in class
x = linspace(0,1,100);
a = sum(Data1);
b = sum(~Data1);
beta1 = betapdf(x,a,b);
figure;
plot(x,beta1)
hold on
% find the mode
[~,idx] = max(beta1);
sprintf('most likely rate is %.3f',x(idx))
% add data from the second experiment
Data2 = [0,0,0,0,0,0,0,1,1,1,1,1]; % data from the second sample in class
a2 = a + sum(Data2);
b2 = b + sum(~Data2);
beta2 = betapdf(x,a2,b2);
plot(x,beta2,'LineWidth',2)
hold on

% find the mode
[~,idx] = max(beta2);
sprintf('most likely rate now becomes %.3f',x(idx))

% if you don't belief in conjugation...let's check it.
prior = beta1;
likelihood = x.^sum(Data2) .* (1-x).^sum(~Data2);
posterior = prior.*likelihood;
[~,idx] = max(posterior);
sprintf('manually get most likely rate: %.3f',x(idx))

hold on
magnifier = max(beta2)/max(posterior);
plot(x,posterior*magnifier,'--','LineWidth',2)

%% FYI: test that the conjugation prior of normal distribution is itself