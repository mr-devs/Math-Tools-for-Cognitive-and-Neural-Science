%% HW # 3 - Question 5 - Math Tools - Matthew DeVerna

clear
close all

%% Question 5: Sampling and Aliasing.

% Load files

load('myMeasurements.mat')

% This file contains:
%   - sig = a vector (length = 120) containing voltage values 
% measure from an EEG electrode (sampled at 100 Hz)

%% Plot "sig" as a function of time, in seconds

N=120                                                           ;
fs = 100                                                        ; % 
buff = .02                                                      ; % An axis buffer so that these plots look a bit cleaner
scaled_time = time/fs                                           ; % The x-values for time
fft_domain = [-N/2: (N/2 -1)] * (fs/N)                          ;
ax_Lims = [0 max(scaled_time)+buff min(sig)-buff max(sig)+buff] ; % Create the axis limits

figure ; % Plot
plot(scaled_time, sig, 'ko-', 'LineWidth', 1.25)
title('EEG Electrode Voltage Over Time')
xlabel({'Time' '(Seconds)'})
ylabel('Voltage')
box off
axis(ax_Lims)

%% (A) Creating and Plotting a Subsample

SubSamplingSummary = {'We can see in the plot below that taking a subsample of a signal is a'
                        'problematic approach to summarize a signal. This type of approach may' 
                        'be linear, because you take an accurate sample at each (equally spaced)'
                        'step, however, this approach will not be shift invariant. We see it is'
                        'not shift invariant because if you sample at every 2nd vs. 3rd vs kth'
                        'element, there is always a potential that you will get a very misleading'
                        'signal.'}

% Plotting a Subsample of Sig at Every Fourth Element

step = 4                                                ; % Every fourth element
step_one = zeros(1,step)                                ; % Create zeros to begin logical array
step_one(4) = 1                                         ; % Fill the fourth element with a one as a selector
steps_in_sig = 120/step                                 ; % We will need this many elements
sub_sample_index = repmat(step_one,[1 ,steps_in_sig])   ; % Create the full array
sub_sample_index = logical(sub_sample_index)            ; % Convert it to a logical array
sig_subsample = sig(sub_sample_index)                   ; % Take the subsample
time_subsample = scaled_time(sub_sample_index)          ; % Take the subsample

figure ; % Plot
plot(scaled_time, sig, 'ko-', 'LineWidth', 1.25)
hold on
plot(time_subsample,sig_subsample, 'r*-', 'LineWidth', 1.25)
hold on
title('EEG Electrode Voltage Over Time')
xlabel({'Time' '(Seconds)'})
ylabel('Voltage')
box off
axis(ax_Lims)
legend('Original Signal ("Sig")','Subsampled Signal', 'Location', 'northeastoutside')

%% Subsampling at Different Rates

% Now we create a subplot to show the difference between sampling at different
% times...

figure ; 
subplot(3,1,1)
plot(scaled_time, sig, 'ko-', 'LineWidth', 1.25)
title('EEG Electrode Voltage Over Time: Full Signal')
ylabel('Voltage')
box off
axis(ax_Lims)

for ii = [2, 3]
step = ii                                               ; % Every fourth element
step_one = zeros(1,step)                                ; % Create zeros to begin logical array
step_one(ii) = 1                                        ; % Fill the fourth element with a one as a selector
steps_in_sig = 120/step                                 ; % We will need this many elements
sub_sample_index = repmat(step_one,[1 ,steps_in_sig])   ; % Create the full array
sub_sample_index = logical(sub_sample_index)            ; % Convert it to a logical array
sig_subsample = sig(sub_sample_index)                   ; % Take the subsample
time_subsample = scaled_time(sub_sample_index)                 ; % Take the subsample

% Plot
subplot(3,1,ii)
plot(time_subsample,sig_subsample, 'bx-', 'LineWidth', 1.25)
title('EEG Electrode Voltage Over Time: Sampling Every ' + string(ii) + ' Elements')
ylabel('Voltage')
box off
axis(ax_Lims)
end

xlabel({'Time' '(Seconds)'})

Plot = {'The bottom two signals are subsamples of the top full signal. You can'
            'clearly see that each one of these subsamples paints a very different'
            'picture as to what the original signal consists of.'}



%% (B) Examine the EEG Results in the Frequency Domain

eeg_amps = abs(fft(sig))            ; % Find the amplitudes of sig
eeg_centered = fftshift(eeg_amps)   ; % Center them around the DC

% Plot the frequency domain / amplitudes
figure ;
stem(fft_domain, eeg_centered, 'filled', 'k-.')
box off
title('EEG Amplitudes in the Frequency Domain')
xlabel({'Frequency' '(Hz)'})
ylabel('Amplitude')

% Find the band with the highest amplitude
max_amp = max(eeg_centered) ; 
max_band = fft_domain(eeg_centered == max_amp)  

LargestBand = {'As we can see the largest response is found in the "beta" band at 30Hz'}

%% (C) Plot Fourier Magnitudes for Downsampled Versions of the Original EEG Signal

FoldingAnswer = {'In the below plots, what we''re seeing is the effect of "folding". Any'
                'frequencies that are higher than the Nyquist frequency in a continuous'
                'signal will be aliased or "folded" into lower frequencies when when we'
                'sample it discretely. One way to handle this is to filter the signal to'
                'remove higher frequencies before the signal is sampled.'}


% Below I plot the original signal as well as each downsampled version to
% show visually the effects of folding within the frequency domain. 
            
figure ; 
subplot(4,1,1)
stem(fft_domain, eeg_centered, 'filled', 'r-.')
box off
title('EEG Amplitudes in the Frequency Domain: Full Sample')
ylabel('Amplitude')

for ii = [2, 3, 4]
step = ii                                                   ; % Every fourth element
step_one = zeros(1,step)                                    ; % Create zeros to begin logical array
step_one(ii) = 1                                            ; % Fill the fourth element with a one as a selector
steps_in_sig = 120/step                                     ; % We will need this many elements
sub_sample_index = repmat(step_one,[1 ,steps_in_sig])       ; % Create the full array
sub_sample_index = logical(sub_sample_index)                ; % Convert it to a logical array
sig_down_sample = zeros(1, length(sig))                     ; % Create array of zeros to fill with samples
sig_down_sample(sub_sample_index) =  sig(sub_sample_index)  ; % Fill this array with the sampled values, creating the downsample
down_sampled_amps = abs(fft(sig_down_sample))               ; % Take the fourier transformed amplitudes
down_sampled_amps_centered = fftshift(down_sampled_amps)    ; % Center them around the DC/zero

% Plot
subplot(4,1,ii)
stem(fft_domain, down_sampled_amps_centered, 'filled', 'b-.')
box off
title('EEG Amplitudes in the Frequency Domain: Downsampled by Every ' + string(ii) + ' Elements')
ylabel('Amplitude')
end
xlabel({'Frequency'})

