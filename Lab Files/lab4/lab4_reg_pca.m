clear
close all

%% regression: geometric intuition in code

% 1. create test data -- for showing the geometry, use only 3 data points
xLuminance = [1:3]'; %Arbitrary units. Column vector of IV
% what is the spike response to different luminance? linear amplification
% plus baseline and noise.
baseline = 2;
noisiness = 3;
amp = 3;
ySpikes = amp * xLuminance + baseline + randn(length(xLuminance),1).*noisiness; %Noisy integrate and fire

%% 2. Plotting this to see what is going on
% 2.1 Plot y against x
figure;
h = plot(xLuminance,ySpikes,'o'); %Put a handle on things
h.MarkerSize = 18;
h.MarkerFaceColor = 'b';
h.MarkerEdgeColor = 'k';
xlabel('Luminance')
ylabel('Spikes')
set(gca,'Fontsize',18)
xlim([0.5 3.5])

%% 2.2 Plot the entire dataset (both IVs and DVs) as vectors
figure
h1 = plot3([0, xLuminance(1)], [0, xLuminance(2)], [0, xLuminance(3)]);
h1.Color = 'b';
h1.LineWidth = 2; 
hold on
h2 = plot3([0, ySpikes(1)], [0, ySpikes(2)], [0, ySpikes(3)]);
h2.Color = 'g';
h2.LineWidth = 2;
legend('Luminance','Spikes')
xlabel('trial 1')
ylabel('trial 2')
zlabel('trial 3')
set(gca,'fontsize',18)
rotate3d on
grid on

%% 3 optimization
% 3.1 vanilla linear regression
% Now let's find the scalar, that - when multiplied with the IV vector,
% minimizes the distance to the measurements. 
% aka we look for beta_op which minimizeds ||y - beta_op * x||^2

% how to find beta_op?
beta_op = ySpikes'*xLuminance / (xLuminance'*xLuminance); %Normalized projection

% with beta_op, how to get the prediction and prediction error?
prediction = beta_op*xLuminance; %Simplest possible case
prederr = ySpikes - prediction;
prederr_op = prederr'*prederr;
% check if the error vector lies perpendicular to the X vector
if round(prederr'*xLuminance,4)==0 % do rounding strictly
    disp('yes, error is orthogonal to the X vector')
end

% Add this to the plot
h3 = plot3 ([0, prediction(1)], [0 prediction(2)], [0 prediction(3)]);
h3.LineWidth = 2;
h3.Color = 'k';
h3.LineStyle = ':';
legend('Luminance','Spikes','Prediction')

% Let's explicitly add the distance between the error
h4 = plot3([ySpikes(1), prediction(1)],[ySpikes(2), prediction(2)], [ySpikes(3), prediction(3)]);
h4.LineWidth = 2;
h4.Color = 'r';
legend('Luminance','Spikes','Prediction','Error')
axis equal

% Convince yourself: is this indeed the best beta?
numBeta = 200; %Fine resolution
testBetas = linspace(0, 10, numBeta);
sqrError = nan*ones(1,numBeta); %initialize to an improbable number (nan = not a number)
for ii = 1:numBeta %Go through all the betas
    prediction = testBetas(ii) .* xLuminance;
    %We need to now explicitly represent the distance metric
    sqrError(ii) = sum((prediction-ySpikes).^2); %Sum of squares
end

% or, non-loop version
predmat = xLuminance*testBetas;
prederrmat = predmat - ySpikes*ones(1,size(predmat,2)); % using vector of 1s to "replicate the vector"
errvec = sum(prederrmat.^2,1);
if sum(round(errvec - sqrError,5))==0
   disp('same error as in for loop')
end

% now check if the prederr_op is indeed optimal?
if min(errvec) >= prederr_op
   disp('seems like we did find a great solution')
end



%Plot that
figure
plot(testBetas,sqrError,'b','linewidth',3)
line([beta_op beta_op], [min(ylim) max(ylim)], 'color','r','linewidth',2);
xlabel('Beta')
ylabel('Error metric')
set(gca,'fontsize',18)
set(gcf,'color','w')

%% 3.2 is this a good model? -- let's add intercept.
% Now we are doing multiple regression!

% So the regression model is Y = beta0 * 1 + beta1*X
% how to write the regression for the problem now?

% Since we generated the data, we know this is close to ground truth.
% baseline and amp parameters should be recovered!
XX = [ones(size(xLuminance)),xLuminance];
[U,S,V] = svd(XX);
disp(diag(S))

% In theory, we wanna minimize err=|Y-X*b|=|U'Y-SV'b|=|Y_star -
% S*b_star|=|Y_star-b_starstar|
% What is the dimensionality for Y_star and b_sstarr?
Y_star = U'*ySpikes;
b_ss = [Y_star(1:length(diag(S)));0];
b_star = b_ss(1:length(diag(S)))./diag(S);
b = V*b_star;
% compare to the original parameters we used to create the spike data...
if sum(round(b - [baseline;amp]))==0
    disp('parameter recovered successfully')
end
disp(b - [baseline;amp]) % why are they not the same sometimes?
% compare the fitting
figure;
plot(xLuminance,XX*b,'b-')
hold on
plot(xLuminance,XX*[baseline;amp],'r-')
legend('fit','ground truth')
hold on
plot(xLuminance,ySpikes,'o','MarkerFaceColor','b') % add raw data

%% 4. go to bigger dataset?
% 4.1 do the fit using the same method; compare with matlab
% let's first write what we have into a function called "lin_reg". 
ndata = 1000;
xLuminance2 = rand(1000,1)*5; % x range in 0 to 5
ySpikes2 = amp * xLuminance2 + baseline + randn(length(xLuminance2),1).*noisiness; %Noisy integrate and fire

% now, get the regressed parameters.
XX2 = [ones(size(xLuminance2)) xLuminance2];
b2 = lin_reg(XX2,ySpikes2);
sqrErrsum = sum((ySpikes2 - XX2*b2).^2);
sprintf('Error = %d',sqrErrsum)
disp (b2 - [baseline;amp])
figure
plot(xLuminance2,XX2*b2,'b-')
hold on
plot(xLuminance2,XX2*[baseline;amp],'r-')
legend('fit','ground truth')
hold on
plot(xLuminance2,ySpikes2,'.') % add raw data

% check with matlab's own function
mdl = fitlm(xLuminance2,ySpikes2);

%% 4.2 landscape of error
% first, create a meshgrid of beta values on the two dimensions
ngrid = 50;
beta0 = linspace(0,5,ngrid);
beta1 = linspace(0,5,ngrid);
% now compute the errors
allerr = nan*ones(ngrid);
for k0 = 1:ngrid
    for k1 = 1:ngrid
        bb = [beta0(k0);beta1(k1)];
        allerr(k0,k1) = sum((XX2*bb-ySpikes2).^2);
    end
end
% check the minimum error should be included
disp([min(allerr(:)),sqrErrsum])
% plot the contour!
figure;
[betaX,betaY] = meshgrid(beta0,beta1);
contour3(betaX,betaY,allerr')
xlabel('beta0')
ylabel('beta1')
zlabel('squared error')
rotate3d on

%{ 
% plot on 2d plane as Mike did on the blackboard
figure;
[betaX,betaY] = meshgrid(beta0,beta1);
contour(betaX,betaY,allerr') % you may also try "contourf"
xlabel('beta0')
ylabel('beta1')
zlabel('squared error')
%}

% plot in the starstar space?
[U,S,V] = svd(XX2);
bss_op = S*V'*b2;
disp(bss_op(1:2))
% sample the values around optimal(aka the center of circle)
samprange = max(abs(bss_op));
bss0 = linspace(-samprange,samprange,ngrid)+bss_op(1); 
bss1 = linspace(-samprange,samprange,ngrid)+bss_op(2);
% now compute the errors
allerr_st = nan*ones(ngrid);
for k0 = 1:ngrid
    for k1 = 1:ngrid
        bb = [bss0(k0);bss1(k1);zeros(length(ySpikes2)-2,1)];
        allerr_st(k0,k1) = sum((bb-U'*ySpikes2).^2);
    end
end
figure;
[betaX,betaY] = meshgrid(bss0,bss1);
contour3(betaX,betaY,allerr_st')
xlabel('beta0_ss')
ylabel('beta1_ss')
zlabel('squared error')
rotate3d on


%% 5. constrained regression, visualized
% let's say we add the constraint that ||beta||=2, meaning it's on a circle
% plot this constraint in both original and starstar space

% what should the constraint look like? 
betamat = [2*cos(linspace(0,2*pi));2*sin(linspace(0,2*pi))];

% what about the constraint in the beta^starstar space?
bssmat = S*V'*betamat;

% let's plot everything together! 
figure;
subplot 211
plot(betamat(1,:),betamat(2,:),'r--','LineWidth',3)
hold on
[betaX,betaY] = meshgrid(beta0,beta1);
ncontour=70; % make the contour denser
contour(betaX,betaY,allerr',ncontour)
xlabel('beta0')
ylabel('beta1')
axis equal

subplot 212
plot(bssmat(1,:),bssmat(2,:),'r--','LineWidth',3)
hold on
[betaX,betaY] = meshgrid(bss0,bss1);
contour(betaX,betaY,allerr_st',ncontour)
xlabel('beta0_ss')
ylabel('beta1_ss')
axis equal

% Bonus question: actually find the solution!

%% B. total least squared regression
% create fake data
Neuron1 = randn(10,1); 
Neuron2 = randn(10,1);
Neuron3 =(Neuron1 + Neuron2)./2; %correlated data

% 1. what's the conceptual difference with the linear regression we just did?
% linear regression: from Neuron1, try to predict the response of Neuron3
figure;
plot(Neuron1,Neuron3,'o')
hold on
beta_op = Neuron3'*Neuron1 / (Neuron1'*Neuron1); %Normalized projection
plot(Neuron1,beta_op*Neuron1,'b-')
hold on
% plot the error
line([Neuron1,Neuron1]',[beta_op*Neuron1,Neuron3]','Color','r')
hold on

% Total least square regression -- see the board
% When we say minimizing ||Du||...what is D, what is u now?
D = [Neuron1,Neuron3];
uvec = [beta_op;-1];
uvec = uvec/sqrt(uvec'*uvec);
line([0,uvec(1)],[0,uvec(2)],'Color','g','Marker','>')
hold on
% plot the error of total least sqr 
proj = D*uvec;
projvec = proj*uvec';
line([Neuron1,Neuron1-projvec(:,1)]',[Neuron3,Neuron3-projvec(:,2)]','Color','g')
axis equal
% what's the total err ||Du||^2?
totalerr = proj'*proj;

% Caution: this is NOT the optimal u! 

% 2. Find minimum total least squared error
[U,S,V] = svd(D);
disp(S)
uvec_min = V(:,2);
totalerr_min = (D*uvec_min)'*(D*uvec_min);
if totalerr<totalerr_min % just a sanity check
    fprintf('wow sth is wrong...minimum err %.3f is not smaller than %.3f',totalerr_min,totalerr)
end
% plot the regression line
regbeta = -uvec_min(1)/uvec_min(2);
figure;
line([-2,2],[-regbeta,regbeta],'Color',[.8,0,.8])
hold on
%plot(Neuron1,Neuron3,'s')
%hold on
% Note we can re-center the data so that the line goes through the center
plot(Neuron1-mean(Neuron1),Neuron3-mean(Neuron3),'o')

% What if we want maximum total squared error instead?
uvec_max = V(:,1);
totalerr_max = (D*uvec_max)'*(D*uvec_max);
if totalerr>totalerr_max % just a sanity check
    fprintf('wow sth is wrong...maximum err %.3f is not bigger than %.3f',totalerr_max,totalerr)
end

%% C. PCA for dimension reduction
% 1. Generate a bunch of neurons that are the linear combination of neuron 1&2
Neuron1 = sin((1:500)*pi/36)'; 
Neuron2 = sin((1:500)*pi/128+pi/4)';
Ntot = 50;
combfac = rand(2,Ntot);
combfac(1,:) = 10*combfac(1,:);
allNeurons = [Neuron1,Neuron2]*combfac; %check the size!

%% 2. Do SVD on this dataset
[U,S,V] = svd(allNeurons);
% What do you expect to see in S?

% let's go check what's in S!
diag(S)

% therefore, the first two rows of V are the directions of biggest total
% square error.
% Let's project the data onto the first three V components and see what we
% get.
figure
plot3(allNeurons*V(:,1),allNeurons*V(:,2),allNeurons*V(:,3),'o')
rotate3d on
xlabel('proj 1')
ylabel('proj 2')
zlabel('proj 3')
zlim([-10,10])

% these components correspond to the actual source signals
figure
subplot 211
plot(Neuron1)
hold on
plot(-allNeurons*V(:,1)/50) % scale it a bit
subplot 212
plot(Neuron2)
hold on
plot(allNeurons*V(:,2))

%% 3. PCA on the variance matrix of the data (check equivalence with SVD)
NeuronsCen = allNeurons - mean(allNeurons);
VarMat = NeuronsCen'*NeuronsCen;
[Vmat, lambdamat] = eig(VarMat);
eigvs = flipud(diag(lambdamat)); % use flip to put big eigen-values at the beginning
% check the eigen values are equal to what we got from SVD
svals = diag(S);
round(sqrt(eigvs(1:2))- svals(1:2))
% check the principle components are equal
figure;
subplot 121
plot(Vmat(:,length(Vmat)),V(:,1),'+')
xlabel('pc1')
ylabel('V1')
axis square
subplot 122
plot(Vmat(:,length(Vmat)-1),V(:,2),'+')
xlabel('pc2')
ylabel('V2')
axis square

%% 4. design your data...ball, pancake or baguette?
ndata = 500;
% ball: all three components unrelated
row1 = rand(1,ndata);
row2 = rand(1,ndata);
row3 = rand(1,ndata);
datatb = [row1;row2;row3]';
datatb_cen = datatb - mean(datatb(:));
figure;
subplot 121
plot3(row1,row2,row3,'.')
subplot 122
eigvals = eig(datatb_cen'*datatb_cen);
bar(eigvals)

% pancake: one components flat
row1 = rand(1,ndata);
row2 = rand(1,ndata);
row3 = rand(1,ndata)*0.1;
%row3 = rand(1,ndata)*0.1;
datatb = [row1;row2;row3]';
datatb_cen = datatb - mean(datatb(:));
figure;
subplot 121
plot3(row1,row2,row3,'.')
zlim([0,1])
subplot 122
eigvals = eig(datatb_cen'*datatb_cen);
bar(eigvals)


% baguette: components on 1 line
row1 = rand(1,ndata);
row2 = row1*.1;
row3 = rand(1,ndata)*0.1;
datatb = [row1;row2;row3]';
datatb_cen = datatb - mean(datatb(:));
figure;
subplot 121
plot3(row1,row2,row3,'.')
zlim([0,1])
subplot 122
eigvals = eig(datatb_cen'*datatb_cen);
bar(eigvals)


%% 5. let's grab some real data...
grdD = table2array(readtable('gradeshw1.csv'));

% each row is a student, each column is a question
% question: are the scores of some questions highly correlated?
table_cen = grdD - mean(grdD);
[Vmat, lambdamat] = eig(table_cen'*table_cen);
figure;
bar(sqrt(diag(lambdamat))) % one component is pretty big! 
% we could project everyone's grad onto this vector...

% let's see the corresponding projection weights of the first 2 components
figure;
subplot 211
stem(Vmat(end,:))
subplot 212
stem(Vmat(end-1,:))
% seems like question 5 and 6 have biggest variance while question 1 & 2
% doesn't vary that much.

% WARNING for PCA: if different dimensions are not of same magnitude, the
% small magnitude dims will always be neglected!
% to save this, normalize every dimension before pca.
table_norm = grdD - repmat(mean(grdD,1),length(grdD),1);
table_norm = table_norm ./ repmat(std(grdD,1),length(grdD),1);
[Vmat, lambdamat] = eig(table_norm'*table_norm);
figure;
bar(sqrt(diag(lambdamat))) 
% let's see the corresponding projection weights of the first 2 components
% again
figure;
subplot 211
stem(Vmat(end,:))
subplot 212
stem(Vmat(end-1,:))
% question 1 & 2 get more weights now.
% so conclusion: all questions are important!!

