%% Homework2 - Question 2 - Matthew DeVerna

% Purpose: This script was written to tackle question # 2 of Math Tools
% homework # 2.

% Author: Matthew DeVerna

% Date: 10/9/19

%% 2. Polynomial Regression

% Find a least-squares fit of the data with polynomials of order 0 (a constant),
% 1 (a line, parameterized by intercept and and slope), 2, 3, 4, and 5. [Compute
% this using svd and basic linear algebra manipulations that you’ve learned in class!] 

% Load the file regress1.mat into your MATLAB environment.

load('regress1.mat')

scatter(x, y)

% Create a matrix for the polynomials
data = [ones(length(x),1), x, x.^2, x.^3, x.^4, x.^5 ] ;

% Take the SVD because we'll need it later
[U,S,V] = svd(data(:,2));

% In theory, we would like to reduce error as much as possible
% err=|Y-X*b| = |U'Y-SV'b| = |Y_star - S*b_star| = |Y_star-b_starstar|

y_star = U'*y ;

% beta_start_star = S# * b_star so we find S_pound
S_pound = diag_P_INVER(S) ;

% Then plug it into the equation as written above (we transform both B-star
% and Y_star by S# - once we do so, we know that beta_star_opt = S#*b_star
beta_star_opt = S_pound * y_star ;

% Finally, we rotate back to the orignial coordinate system
beta = V*beta_star_opt ;

% All equations below are equivalent

beta_opt_1 = V*beta_star_opt
beta_opt_1a = V * S_pound * y_star
beta_opt_1b = V * S_pound * U' * y

% You can also do the below for only one regressor
the_short_cut_method = y'*x/(x'*x)

sprintf('As you can see, all three forms of the equation return the same value for beta_opt.\n')

% Lets create a function for plotting later
beta_opt_Func1 = x*beta_opt_1 ;

% And repeat, using hte linear_reg() function for all other polynomials

line = data(:,1:2) ;
beta_opt2 = linear_Reg(line,y)

order2 = data(:,1:3) ;
beta_opt3 = linear_Reg(order2 ,y)

order3 = data(:,1:4) ;
beta_opt4 = linear_Reg(order3,y)

order4 = data(:,1:5) ;
beta_opt5 = linear_Reg(order4,y)

order5 = data(:,1:6) ;
beta_opt6 = linear_Reg(order5,y) 

% Create functions for all of these as well...

beta_opt_Func2 = line*beta_opt2     ;
beta_opt_Func3 = order2*beta_opt3   ;
beta_opt_Func4 = order3*beta_opt4   ;
beta_opt_Func5 = order4*beta_opt5   ;
beta_opt_Func6 = order5*beta_opt6   ;

% Plot the data and fit...

% Just the regressors
subplot(2,3,1)
scatter(x,y)
hold on
plot(x, beta_opt_Func1, 'b')
title('Constant')

% Regressors plus constant, creating a line
subplot(2,3,2)
scatter(x,y)
hold on
plot(x, beta_opt_Func2, 'r')
title('Line')

% Second Order
subplot(2,3,3)
scatter(x,y)
hold on
plot(x, beta_opt_Func3, 'g')
title('Second Order')

% Third Order
subplot(2,3,4)
scatter(x,y)
hold on
plot(x, beta_opt_Func4, 'c')
title('Third Order')

% Fourth Order
subplot(2,3,5)
scatter(x,y)
hold on
plot(x, beta_opt_Func5, 'k')
title('Fourth Order')

% Fifth Order
subplot(2,3,6)
scatter(x,y)
hold on
plot(x, beta_opt_Func6, 'k')
title('Fifth Order')

% Create a legend for each subplot
for ii = 1:6
    subplot(2,3,ii)
    legend('Data', 'Fit', 'Location', 'northwest')
end

%% Finding Squared Error and Plotting that

% Below we are taking the total squared difference and setting it to a
% specific variable for later plotting...

constant_SE = sum((y - beta_opt_Func1).^2)
line_SE = sum((y - beta_opt_Func2).^2)
order2_SE = sum((y - beta_opt_Func3).^2)
order3_SE = sum((y - beta_opt_Func4).^2)
order4_SE = sum((y - beta_opt_Func5).^2)
order5_SE = sum((y - beta_opt_Func6).^2)

% Create a matrix of total squared error
all_sqrd_err = [constant_SE, line_SE, order2_SE, order3_SE, order4_SE, order5_SE] ;

% Create labels for the below bar plot
orders = {'Zero Order' 'First Order' 'Second Order' 'Third Order' 'Fourth Order' 'Fifth Order'} ;
% Re-order them so they are decrease in error
orders = reordercats(categorical(orders'), {'Zero Order' 'First Order' 'Second Order' 'Third Order' 'Fourth Order' 'Fifth Order'}) ;

% Open figure and plot
figure
bar(orders, all_sqrd_err)
title('Squared Error as a Function of Polynomial Fit')
ylabel('Squared Error')
xlabel('Polynomial')
%% Best Fit Answer

answer = sprintf('Reviewing the plot of squared error we can see that the largest reduction in error comes as we move from the first order to the second order\npolynomial. Additionally, we see a relatively large reduction in error with the third order polynomial as well.\n However, we can see at the end of this fit line in the negative space (on the left), that the line begins to dip back down as the data continues upward.\n\n Noticing this, I would choose the second order fit. It captures the general shape/trend of the data without the need for another parameter which may not generalize.')














