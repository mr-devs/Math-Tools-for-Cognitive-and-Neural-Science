%% Homework2 - Question 4 - Matthew DeVerna

% Purpose: This script was written to tackle question # 2 of Math Tools
% homework # 4.

% Author: Matthew DeVerna

% Date: 10/9/19

%%  Load File

load('windowedSpikes.mat')

% data = 400x150 matrix of whos rows are electrode measurements ( voltages
% recorded for each 150 msec window, at 1msec intervals). 

%% (A) Plot the waveforms superimposed

figure
plot(data')
xlabel('Time')
ylabel('Voltage')
legend('Recordings')
title('Dr. Bell and Dr. Zell''s Electrode Measurement Data')

% I cannot think of a scientifically valid procedure to determine how many
% neurons produced these spikes (only using the plot). We must used PCA/SVD to find the
% eigenvectors/eigenvalues and determine which dimensions are responsible
% for the most variance


%% (B) Plot the Eigen Values in Descending Order

% We can use the svd to find the eigen values (diag(S).^2) 
% and eigen vectors (V')

% Compute the SVD
[U,S,V] = svd(data)     ;

% Since we know that our eigen values are the diagonal components of S,
% squared. We can take these directly from S
eig_vals = (diag(S)).^2 ;

% Now we plot these values in a simple bar graph and it's quite clear
semilogy(eig_vals, 'o') % Make y-axis log
xlabel('Components')
ylabel('Values')
title('Eigen Values for Electrode Data')
legend('Eigen Values')
grid on
box off


%% (C) Project onto the top two princple components

% We know that V' represents our eigen vectors so we take these and utilize
% them to project our data these dimensions

eigen_vectors = V'          ;

% Find our principle components
PC_1 = eigen_vectors(1,:)   ;
PC_2 = eigen_vectors(2,:)   ;
PC_3 = eigen_vectors(3,:)   ; % Find the third for later use

% Project our data onto the new coordinate system
X = data * PC_1'            ;
Y = data * PC_2'            ;
Z = data * PC_3'            ; % Find third for later
x_zeros = zeros(400,1)      ; % Create zeros for 3D vector plotting.

% Plot in the 2D space
plot(X,Y, 'x')
title('Data Projected on Two Largest Principle Components')
xlabel('Principle Component #1')
ylabel('Principle Component #2')
legend('Recordings')
grid on 
box off

%% (D) Project onto the top three principle components

plot3(X,Y,Z, 'x')
rotate3d on 
title('Data Projected on Three Largest Principle Components')
xlabel('Principle Component #1')
ylabel('Principle Component #2')
zlabel('Principle Component #3')
legend('Recordings')
grid on 
box off
