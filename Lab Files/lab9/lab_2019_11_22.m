%% Bootstrap: Aspirin
clear
aspirin_heart = 104;
aspirin_total = 11037;

placebo_heart = 189;
placebo_total = 11034;

% like flipping a coin...
asp_data = zeros(aspirin_total,1); % 11037 total subjects, all zero
asp_data(1:aspirin_heart) = 1; % set 104 subjects to 1 (1 = heart attack)

placebo_data = zeros(placebo_total, 1); % 11034 tot subjects, all zero
placebo_data(1:placebo_heart) = 1; % 189 had heart attacks

ratio_empirical = (aspirin_heart/aspirin_total)/(placebo_heart/placebo_total);
n_boot = 10000;
ratio_boot = zeros(n_boot, 1);
for boot = 1:n_boot
    boot_asp = randsample(asp_data, aspirin_total, 'true');
    boot_placebo = randsample(placebo_data, placebo_total, 'true');    
    n_boot_asp = sum(boot_asp);
    n_boot_placebo = sum(boot_placebo);    
    ratio_boot(boot) = (n_boot_asp/aspirin_total)/(n_boot_placebo/placebo_total);        
end

figure
close all
hold on
histogram(ratio_boot, 'edgecolor', 'none')
xlabel('ratio')
ylabel('count')
plot([1, 1]*ratio_empirical, [0, 1]*700,'--')

%% Bootstrapped model fitting:
clear
rng(2019)
n_pts = 20;
X = (1:n_pts)';
Y = (1*X) + 2 + 3*randn(n_pts,1); % Y= mx +b + noise

close all
figure
subplot(131)
scatter(1:n_pts, Y,'k')
axis square
subplot(132)
scatter(1:n_pts, Y,'k')
hold on

model = [ones(n_pts,1) X];
n_boot = 500;
boot_fit = zeros(n_boot, n_pts); % create empty matrix for fits
boot_beta = zeros(n_boot,2); % create empty matrix for fit coeffs

for boot=1:n_boot

boot_ind = randi(n_pts, n_pts, 1); % index to resample data
boot_model = model(boot_ind,:); % subsample rows of model
boot_Y = Y(boot_ind); % subsample % subsample corresponding Y- data

beta = regress(boot_Y, boot_model); % fit this boostrapped data
boot_beta(boot,:) = beta; % store the betas

boot_fit(boot,:) = model*beta; % get the fit by applying beta to model
plot(boot_fit(boot,:),'b')

end

axis square

% compute the standard deviation at each point
boot_fit_std = std(boot_fit,[],1);

% empirical non-bootstrapped fit
beta_emp = regress(Y,model);

emp_fit = model*beta_emp;

subplot(133)
hold on
scatter(1:n_pts, Y ,'k')
errorbar(1:n_pts, emp_fit, boot_fit_std ,'capsize',0)
axis square

figure
subplot(121)
histogram(boot_beta(:,2))
title('slope')
subplot(122)
histogram(boot_beta(:,1))
title('intercept')

%% Permutation test: paired t-test
clear
close all
n_subj = 20;
drug_before = 100 + 5*randn(n_subj,1); % IQ
drug_after = drug_before + 3.2 + 2*randn(n_subj,1); % post-drug

figure
plot([drug_before, drug_after]' ,'k' ,'linewidth', 2)
set(gca,'xtick',1:2, 'xticklabel',{'before','after'})
xlim([0.9,2.1])
ylabel('IQ')

emp_mean_diff = mean(drug_after-drug_before);

n_perm = 10000;
all_data = [drug_before, drug_after]; % concatenate the data
perm_diffs = zeros(n_perm,1);
for perm = 1:n_perm
   perm_ind = randperm(2*n_subj); % generate random indices
   
   perm_before= all_data(perm_ind(1:n_subj))';
   perm_after = all_data(perm_ind(n_subj+1:end))';
   
   perm_diffs(perm) = mean(perm_after-perm_before);      
    
end

figure
plot([perm_before perm_after]','k' ,'linewidth', 2)
set(gca,'xtick',1:2, 'xticklabel',{'before','after'})
xlim([0.9,2.1])
ylabel('IQ')

figure
hold on
histogram(perm_diffs,'edgecolor','none')
plot([1 1]*emp_mean_diff, [0,500],'r--','linewidth',3)
xlabel('Mean difference')
ylabel('Count')
legend('Null', 'Empirical')

% Compute p value
% how many null values lie BELOW the real value? 
% Is the empirical value greater than 95% of the null?

p_val = 1 - sum(emp_mean_diff>perm_diffs)/n_perm;

%% Permutation test: continuous data fitting
clear
close all
rng(2019)
n_pts = 20;
X = (1:n_pts)';
Y = (1*X) + 2 + 3*randn(n_pts,1); % Y= mx +b + noise

model = [ones(n_pts,1) X];
beta = regress(Y,model);

figure
subplot(121)
hold on
cmap = parula(n_pts);
scatter(X,Y,50,cmap,'filled')
plot(X,model*beta)
legend('data','fit')
axis square
disp(beta) 

% are these beta fit coeffs different than chance?

n_perm = 10000;
beta_perm = zeros(n_perm,2);
for perm = 1:n_perm
    perm_ind = randperm(n_pts);
    X_perm = X(perm_ind);
    
    model_perm = [ones(n_pts,1) X_perm];
    beta_perm(perm,:) = regress(Y, model_perm);    
    
end

subplot(122)
hold on
scatter(X_perm,Y,50,cmap,'filled')
plot(X_perm, model_perm*beta_perm(end,:)')
axis square

figure
subplot(121)
hold on
histogram(beta_perm(:,2),'edgecolor','none') % slope
plot([1 1] * beta(2), [0 1]*500,'--','linewidth',3)
legend('10k shuff null', 'empirical')
ylim([0,800])
title('Slope')

subplot(122)
hold on
histogram(beta_perm(:,1), 'edgecolor','none')  % intercept
plot([1 1] * beta(1), [0 1]*500,'--','linewidth',3)
ylim([0,800])
title('Intercept')

% How do we compute p-values?
% Find how many values of shuffled slope coeffs are higher/lower than true
% value

%% Permutation test: EEG data spectrogram
clear all
close all
load eeg.mat

imagesc(eeg)
set(gca,'xtick',[],'ytick',[])
xlabel('Time')
ylabel('Frequency')
h=colorbar;
axis square 
ylabel(h,'Relative power')


n_perm = 5000;
perm_eeg = zeros(size(eeg,1),size(eeg,2),n_perm); % Create null matrix
n_eeg = numel(eeg);
for perm = 1:n_perm    
    perm_ind = randperm(n_eeg,n_eeg);
    perm_eeg(:,:,perm) = reshape(eeg(perm_ind),size(eeg));    
end

% Plot three examples
figure
for plt = 1:3
    subplot(1,3,plt)   
    imagesc(perm_eeg(:,:,plt))
    set(gca,'ytick',[],'xtick',[])
    axis square
end

figure
hold on
histogram(perm_eeg(1,1,:)) % grab every null value at the first pixel
plot([1, 1]*eeg(1), [0,1] * 600,'--','linewidth',2) % plot the actual value
legend('Null distribution for pixel 1','Empirical value')

% Find 2.5 and 97.5 percentile for each pixel
thresh_low = prctile(perm_eeg, 2.5, 3);  % find min along 3rd dimensin
thresh_high = prctile(perm_eeg, 97.5, 3); % find max

subplot(121)
imagesc(thresh_low,[min(eeg(:)), max(eeg(:))])
colorbar
title('lower bound')
axis square
set(gca,'xtick',[],'ytick',[])

subplot(122)
imagesc(thresh_high, [min(eeg(:)), max(eeg(:))])
colorbar
title('upper bound')
axis square
set(gca,'xtick',[],'ytick',[])

% plot the thresholded image
% These are the values that are statistically different than the null
mask = eeg<thresh_low | eeg>thresh_high;
thresholded = eeg .* mask;
imagesc(thresholded, [min(eeg(:)), max(eeg(:))])

figure
colorbar
axis square
set(gca,'xtick',[],'ytick',[])
xlabel('Time')
ylabel('Frequency')
ylabel(h,'Relative power')

%% Cross-validation -- model fitting
clear all
load regress1.mat
n_pts = size(y,1);
XX = [ones(n_pts,1), x, x.^2,x.^3,x.^4,x.^5];
n_models = size(XX,2);

split_ind = randperm(n_pts);

figure
scatter(x,y)

hold on
n_folds = n_pts; 

mse_xval = zeros(n_folds, n_models);

% There are lots of ways to do x-validation, but leave-one-out is simplest
% to implement & understand

for mdl = 1:n_models
    for fold = 1:n_pts % Leave-one-out cross-validation
        
        xval_ind = true(n_pts,1);
        xval_ind(fold)=false;
        
        x_train = XX(xval_ind,1:mdl);
        y_train = y(xval_ind); 
        
        x_test = XX(~xval_ind,1:mdl);
        y_test = y(~xval_ind);
        
        % compute the fit
        [U,~,V] = svd(x_train);
        S_vec = svd(x_train); %get vector of singular values
        S_inv = zeros(size(x_train));
        S_inv(1:mdl,1:mdl) = diag(1./S_vec);
        b = V*S_inv'*U'*y_train; 
        
        fit_train =  XX(:,1:mdl)*b; % trained model
        fit_test = x_test * b; % model's prediction for test point
        
        mse_xval(fold,mdl) = mean((y_test - fit_test).^2); 
                
        plot(x, fit_train)
    end
end

% Plot the MSE and error bars 
% We want the model that minimizes error AND is the simplest
% The simplest model that lie within 2 standard errors of the minimum MSE
% will be chosen as the 'best' one.

figure
errorbar(0:5, mean(mse_xval), 2*sqrt(var(mse_xval)/n_pts),'capsize',0,'linewidth',3)
xlabel('Polynomial order')
ylabel('MSE')

hold on
% Find LOWEST MSE
[min_mse, min_ind] = min(mean(mse_xval)); % lowest MSE is 5th index, i.e. 4th order polynomial

% Then find 95% confidence interval upper threshold of that MSE, the SIMPLEST model
% whose MSE confidence interval lies within this threshold is the model we
% will go with
thresh_mse = min_mse + 2*sqrt(var(mse_xval(:,min_ind)))/n_pts;

plot(0:5, ones(n_models,1)* thresh_mse ,'k--')
legend('MSE and 95% conf. int.', 'best fit thresh')
set(gca,'xtick',0:5)
% Get the index of the smallest order polynomial that lies within
% thresh_mse

lower_bound = mean(mse_xval) -  2*sqrt(var(mse_xval)/n_pts);

lower_bound<thresh_mse

disp('Models polynomial order 3, 4, 5 all lie within the threshold')
disp('Model w/ poly order 3 is the simplest, so that is the best cross-validated model')

figure
hold on
scatter(x,y)
beta_final = regress(y,XX(:,1:4)); % we found taht 3rd order poly is best
plot(x,XX(:,1:4)*beta_final)
legend('Data', 'Best fit model, assessed using LOO cross-validation')
