%% MATH TOOLS LAB 5: CONVOLUTION

%% 1 Modeling event-related fMRI BOLD signal

% Load the HRF for convolution
load('hrf.mat'); 

% Plot the hrf to see what it looks like
figure
plot(hrf)
title('Canonical HRF')
xlabel('time (s)')
ylabel('BOLD')

% What do you notice about the signal?

% Quick aside: sampling at the standard 2s TR
figure
hold on
plot(hrf)
scatter(2:2:19, hrf(2:2:19))
title('Canonical HRF, TR = 2s')
xlabel('time (s)')
ylabel('BOLD')

% Assume we have a very brief stimulus that we predict will elicit 
% activity in the voxel we care about. We can represent this activity with
% a single 1 in a vector of zeros.

taskDesign = zeros(1,100);
eventTR = 30;
taskDesign(eventTR) = 1;

% The neural activity will cause the relative levels of deoxygenated blood
% in the voxel to change, which is reflected in the BOLD response. 

predResp = convn(taskDesign, hrf, 'same');  % predicted response

figure
hold on
plot(predResp)
stim = plot([eventTR, eventTR], [0 1], 'k');
title('Impulse and corresponding response')
xlabel('time (s)')
ylabel('BOLD')
legend(stim,{'stimulus'})

% That doesn't look like an event-triggered response...why? Where does 
% the HRF appear in the predicted response?

delay = eventTR - ceil(length(hrf)/2);

figure
for ii = 1:length(hrf)
    
    % Plot kernel and incrementally highlight each point
    subplot(1,2,1)
    plot(hrf)
    hold on
    title('Canonical HRF')
    kernel_point = scatter(ii,hrf(ii),'r');
    xlabel('time (s)')
    ylabel('BOLD')
    legend(kernel_point, 'kernel')
    hold off
    
    % Plot response and incrementally highlight corresponding point
    subplot(1,2,2)
    plot(predResp)
    hold on
    stim = plot([eventTR, eventTR], [0 1], 'k');
    title('Impulse and corresponding response')
    delayed_kernel = scatter(delay+ii, predResp(delay+ii),'r');
    xlabel('time (s)')
    ylabel('BOLD')
    legend([stim,delayed_kernel],{'stimulus','delayed kernel point'})
    hold off
    shg
end

% Past research has repeatedly demonstrated that the peak should be
% ~6 seconds after the stimulus. How can we change the kernel so that our
% predicted response is actually what we would predict?

paddedhrf = [zeros(1,20) hrf];
newPredResp = convn(taskDesign, paddedhrf, 'same');

% Now plot our predicted response
figure
plot(newPredResp)
hold on
stim = plot([eventTR, eventTR], [0 1], 'k');
peak = scatter(eventTR+6, newPredResp(eventTR+6),'r');  % plot point 6 seconds after onset
title('Impulse and corresponding response (aligned)')
xlabel('time (s)')
ylabel('BOLD')
legend([stim peak],{'stimulus', 'peak'})


%% Exercise: time-varying trials
% Add more trials to the experiment, and see what happens when you vary
% how close or far apart they are in time. Plot the corresponding
% predicted responses and note how they change.

% Here's an example for three trials, each 20 seconds apart.

taskDesign2 = zeros(1,100);
eventTR = 30;
taskDesign2(eventTR) = 1;
eventTR1 = 10;
taskDesign2(eventTR1) = 1;
eventTR2 = 50;
taskDesign2(eventTR2) = 1;

paddedhrf = [zeros(1,20) hrf];
newPredResp = convn(taskDesign2, paddedhrf, 'same');

% Plot our predicted response
figure
plot(newPredResp)
hold on
stim = plot([eventTR, eventTR], [0 1], 'k');
plot([eventTR1, eventTR1], [0 1], 'k')
plot([eventTR2, eventTR2], [0 1], 'k')
peak = scatter(eventTR+6, newPredResp(eventTR+6),'r');
scatter(eventTR1+6, newPredResp(eventTR1+6),'r')
scatter(eventTR2+6, newPredResp(eventTR2+6),'r')
title('Impulse and corresponding response for more trials')
xlabel('time (s)')
ylabel('HRF')
legend([stim peak],{'stimulus', 'peak'})

%% 2 Convolution to denoise LFP data

% Load the LFP data
load('LFPdata.mat');

% Now the workspace will contain a variable "A" that contains
% voltage samples per second, but only for 3000 milliseconds. They
% represent the "local field potential", voltages in microvolt at the
% electrode tip. 

% Plot the data
figure
plot(A)
title('Raw Local Field Potential')
xlabel('time (ms)')
ylabel('voltage')

% By visual inspection, we note the presence of a signal, but also of a lot
% of noise. The sharp stuff in particular does not look very
% physiological. There are no known processes in the brain that could
% generate such signals, but there are a lot of known problems with
% electrical noise that could. Solution: Smoothing by convolution.


%% Exercise: kernel size

% How big should we make the kernel to minimize noise and keep our signal? 
% Try a bunch of kernel sizes and convolve them with the input to see what
% looks best (hint: use the simpest kernel you can think of that will
% average out noise and increase its size incrementally).

% What happens as the kernel gets larger?
%   - to the noise
%   - to the signal
%   - to the length of the output

for ii = 1:length(A)                % go all the way
    kernelSize = ii;                % this is the length of any given kernel
    kernel = ones(kernelSize,1);    % symmetric, equal weights, the simplest one
    kernel = kernel./sum(kernel);   % normalize so that the output values don't scale with kernel length
    pause(0.15)
    Aconv = conv(A,kernel,'valid');
    plot(Aconv)
    title(['Kernel size = ', num2str(ii), ' samples'])
    ylim([-40 20])                  % keep y axis constant
    xlim([-100 3100])               % keep x axis constant
    shg
end


%% 3 2D image convolution and boundary handling

% Let's create an artifical image in order to really see the effects of 2D
% convolution. We can do this with luminance ramps.

numSteps = 7;                               % num luminance steps
minLum = 0;                                 % lowest luminance value
maxLum = 255;                               % highest luminance value in an 8 bit system
temp = linspace(minLum, maxLum, numSteps);  % determine the step size
stepSize = unique(round(diff(temp),4));

% To create an artificial 2D image, we can use the funciton meshgrid. 
% Its 2 outputs are x and y, which are a full cross of the inputs, keeping
% cone dimension constant while incrementing the other and vice versa. 
[x,y] = meshgrid(minLum:stepSize:maxLum,minLum:stepSize:maxLum);

% Now plot these to visualize
figure
subplot(1,2,1)
imagesc(x)
subplot(1,2,2)
imagesc(y)
colorbar
colormap(bone)

% Matlab represents images as 3D matrices. The first 2 dimensions are x and
% y coordinates, the third dimension are R, G and B sheets. The values
% represent luminances in 8 bit, so they range from 0 to 255.

clear IM                              % this will be a 3D matrix that represents the image. We create the sheets one by one
IM(:,:,1) = uint8(x);                 % red channel (x)
IM(:,:,2) = uint8(zeros(length(x)));  % silence the green channel to get purple
IM(:,:,3) = uint8(y);                 % uint8 interprets/renders the values in 8 bit, as an unsigned integer.

% plot our image that gets decreasingly red in the x-dimension and
% decreasingly blue in the y-dimension
figure
image(IM)
axis square
title('Unfiltered')


%% Exercise: convolution and boundary handling
% Create a 4x4 averaging kernel and convolve it with the original image. 
% You should be able to clearly see that the convolved image has less
% squares, which are each an average of the colors around them. Then,  
% repeat the convolution with different boundary handling conditions. 
% Reason about why you get the results that you do.

% Create an averager kernel using valid boundary handling
kernel = ones(4);                                   % small summing kernel 
kernelWeight = sum(kernel(:));              
kernel = kernel ./ kernelWeight;                    % normalize to make it an averager
convIM = uint8(round(convn(IM, kernel, 'valid'))); 

figure
image(convIM)
axis square
title('Convolved (valid)')

% Now repeat for full and same boundary handling
convIM = uint8(round(convn(IM, kernel, 'full'))); 
figure
image(convIM)
axis square
title('Convolved (full)')

convIM = uint8(round(convn(IM, kernel, 'same'))); 
figure
image(convIM)
axis square
title('Convolved (same)')


%% 4 Image processing with natural images

% Read in a test image
img = imread('bike.jpg');
imshow(img)

% Blurring filter: convolve the image with a summing kernel 
figure
subplot(2,1,1)
imshow(img)
title('Original Image')

subplot(2,1,2)
kernelWidth = 10; 
kernel = ones(kernelWidth);
kernelWeight = sum(kernel(:));
kernel = kernel/kernelWeight;               % normalize
blur = uint8(convn(img, kernel, 'same'));
imshow(blur)
title('Blurred Image')

% What do you notice about the blurred image? How would you adjust the
% amount that we blur?

% Sharpening filter: contrast enhancement with a difference kernel
figure
subplot(2,1,1)
imshow(img)
title('Original Image')

% Create a contrast-enhancing filter from scratch. Here's an example:
kernel = [];
kernel(1,:) = [-1/9 -1/9 -1/9]; 
kernel(2,:) = [-1/9  1   -1/9]; 
kernel(3,:) = [-1/9 -1/9 -1/9]; 
kernelWeight = sum(kernel(:)); 
kernel = kernel ./ kernelWeight;            % normalize
subplot(2,1,2)
sharper = uint8(convn(img,kernel,'same'));
imshow(sharper)
title('Sharper Image')

% What do yout notice about the sharpened image? How would you adjust the
% amount that we sharpen?


%% Exercise: Sobel operator

% A famous example of 2D image convolution in computer vision is called 
% a Sobel operator, which is a specific type of filter used for edge 
% detection. It works by computing the approximate gradient at each pixel 
% in the image in both the x and y direction, thereby emphasizing edges. 
% Look up the kernel(s) for this operator, and compute both the
% vertical and horizontal gradients (hint: you don't need to normalize in 
% this case, just apply the filter as is). Visualize each, and then
% combine the two images to see the full set of edges.

% Define the kernel for the x direction (y is just the transpose)
sobel_kernel = [];
sobel_kernel(1,:) = [-1 0 1]; 
sobel_kernel(2,:) = [-2 0 2]; 
sobel_kernel(3,:) = [-1 0 1]; 

% Compute the vertical gradient and plot
figure
subplot(1,3,1)
sobelx = uint8(convn(img,sobel_kernel,'same'));
imshow(sobelx)
title('Sobel X')

% Redo above for the horizontal gradient and plot
subplot(1,3,2)
sobely = uint8(convn(img,sobel_kernel','same'));
imshow(sobely)
title('Sobel Y')

% Put the images together for full edge detection
subplot(1,3,3)
sobel = sobelx + sobely;
imshow(sobel)
title('Sobel')

