%% HW # 3 - Question 1 - Math Tools - Matthew DeVerna

clear all
close all

%% Question 1: LSI Systems Characterization

% We are trying to experimentally characterize three auditory neurons, in
% terms of their responses to sounds. 

% For this problem, the responses are contained in the compiled matlab
% functions, "unknownSystemX.p" (X = 1,2,3). 
% INPUT - a column vector with length N = 64. (Represents sound pressure over
% time).
% OUTPUT - a column vector of the same length. (Represents mean spike count
% over time). 

% The task is examine each neuron and see if it behaves like it is:
    % i) Linear
    % ii) Shift invariant 
    % iii) Handling boundaries in a circular manner
    
%% A) "Kick the Tires"

%% pt. 1 Checking for Shift Invariance
% Measure the response to an impulse in the first position of an input
% vector. Check that this impulse response is shift-invariant by comparing
% to the response an impulse at positions n = 2,4,8.

% Create general impulses shifted by some amount
zeros_impulse = zeros(64,1) ;
impulse = zeros_impulse ;
impulse(1) = 1 ;
impulse_2 = circshift(impulse, 1) ;
impulse_4 = circshift(impulse, 3) ;
impulse_8 = circshift(impulse, 7) ;

% Create System 1 impulse responses at each position
imp_resp_1_1 = unknownSystem1(impulse)    ;
imp_resp_1_2 = unknownSystem1(impulse_2)  ;
imp_resp_1_4 = unknownSystem1(impulse_4)  ;
imp_resp_1_8 = unknownSystem1(impulse_8)  ;

% Check that these shifts are equal by taking their difference
shift_check1_1 = [imp_resp_1_2 - circshift(imp_resp_1_1,1) < 1e-8]'
shift_check1_2 = [imp_resp_1_4 - circshift(imp_resp_1_1,3) < 1e-8]'
shift_check1_3 = [imp_resp_1_4 - circshift(imp_resp_1_2,2) < 1e-8]'
shift_check1_4 = [imp_resp_1_8 - circshift(imp_resp_1_1,7) < 1e-8]'

% Print Answer
System_1_ShiftCheck1 = {'We can see that system 1 is shift invariant because we get the same'
    'response whether we provide an input at position 1, 2, 4 or 8 - the'
    'response is simply shifted by the position of the impulse position.'} 

%% Plot the Responses Visually
system1 = figure                    ;
resp1 = subplot(4,1,1)              ;
stem(imp_resp_1_1, 'b', 'filled')   ;
title('Position 1 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;
resp2 = subplot(4,1,2)              ;
stem(imp_resp_1_2, 'r', 'filled')   ;
title('Position 2 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;
resp3 = subplot(4,1,3)              ;
stem(imp_resp_1_4, 'k', 'filled')   ;
title('Position 4 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;
resp4 = subplot(4,1,4)              ;
stem(imp_resp_1_8, 'g', 'filled')   ;
title('Position 8 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;
sgtitle('System 1 Shift Invariance Check Responses') ;

%% Find System 2 impulse responses at each position
imp_resp_2_1 = unknownSystem2(impulse)    ;
imp_resp_2_2 = unknownSystem2(impulse_2)  ;
imp_resp_2_4 = unknownSystem2(impulse_4)  ;
imp_resp_2_8 = unknownSystem2(impulse_8)  ;

% Check that these shifts are equal by taking their difference
shift_check2_1 = [imp_resp_2_2 - circshift(imp_resp_2_1,1) < 1e-8]'
shift_check2_2 = [imp_resp_2_4 - circshift(imp_resp_2_1,3) < 1e-8]'
shift_check2_3 = [imp_resp_2_4 - circshift(imp_resp_2_2,2) < 1e-8]'
shift_check2_3 = [imp_resp_2_8 - circshift(imp_resp_2_1,7) < 1e-8]'

% Print Answer
System_2_ShiftCheck1 = {'We can see that system 2 is also shift invariant because we get the same'
    'response whether we provide an input at position 1, 2, 4 or 8 - the'
    'response is simply shifted by the position of the impulse position.'} 

%% Plot the responses visually
system2 = figure                    ;
resp1 = subplot(4,1,1)              ;
stem(imp_resp_2_1, 'b', 'filled')   ;
title('Position 1 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;
resp2 = subplot(4,1,2)              ;
stem(imp_resp_2_2, 'r', 'filled')	;
title('Position 2 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;
resp3 = subplot(4,1,3)              ;
stem(imp_resp_2_4, 'k', 'filled')   ;
title('Position 4 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;
resp3 = subplot(4,1,4)              ;
stem(imp_resp_2_8, 'g', 'filled')   ;
title('Position 8 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;
sgtitle('System 2 Shift Invariance Check Responses') ;

%% Find the impulse response at each position for System 3
imp_resp_3_1 = unknownSystem3(impulse)    ;
imp_resp_3_2 = unknownSystem3(impulse_2)  ;
imp_resp_3_4 = unknownSystem3(impulse_4)  ;
imp_resp_3_8 = unknownSystem3(impulse_8)  ;

% Check the system is shift invariant 
shift_check3_1 = [imp_resp_3_2 - circshift(imp_resp_3_1,1) < 1e-8]'
shift_check3_2 = [imp_resp_3_4 - circshift(imp_resp_3_1,3) < 1e-8]'
shift_check3_3 = [imp_resp_3_4 - circshift(imp_resp_3_2,2) < 1e-8]'
shift_check3_4 = [imp_resp_3_8 - circshift(imp_resp_3_1,7) < 1e-8]'

% Print Answer
System_3_ShiftCheck1 = {'We can see that system 3 is NOT shift invariant because the'
    'responses we get are different responses at certain points in the impulse response.'
    'This tells us that the response VARIES as we shift along this axis.'} 

%% Plot responses visually
system3 = figure                    ;
resp1 = subplot(4,1,1)              ;
stem(imp_resp_3_1, 'b', 'filled')	;
title('Position 1 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;       
resp2 = subplot(4,1,2)              ;
stem(imp_resp_3_2, 'r', 'filled')	;
title('Position 2 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;
resp3 = subplot(4,1,3)              ;
stem(imp_resp_3_4, 'k', 'filled')	;
title('Position 3 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;
resp3 = subplot(4,1,4)              ;
stem(imp_resp_3_8, 'g', 'filled')	;
title('Position 3 Responses')       ;
ylabel('Mean Spike Count')          ;
xlabel('Time')                      ;
sgtitle('System 3 Shift Invariance Check Responses') ;

%% pt. 2 Checking for Additivity

% Here we will check whether a system is linear by testing for additivity. We
% will take multiple combinations of initial impulse vectors and see if
% these combinations elicit the same response as the summation of their
% individual outputs, combined.

% First create new combined input vectors via summation
imp_1and2 = impulse + impulse_2     ;
imp_1and4 = impulse + impulse_4     ;
imp_2and4 = impulse_4 + impulse_2   ;
imp_2and8 = impulse_8 + impulse_2   ;

% Then I find the system responses based on these combined 
% input vectors

% System 1
S1_combined_imp_1_2 = unknownSystem1(imp_1and2)    ;
S1_combined_imp_1_4 = unknownSystem1(imp_1and4)    ;
S1_combined_imp_2_4 = unknownSystem1(imp_2and4)    ;
S1_combined_imp_2_8 = unknownSystem1(imp_2and8)    ;
% System 2
S2_combined_imp_1_2 = unknownSystem2(imp_1and2)    ;
S2_combined_imp_1_4 = unknownSystem2(imp_1and4)    ;
S2_combined_imp_2_4 = unknownSystem2(imp_2and4)    ;
S2_combined_imp_2_8 = unknownSystem2(imp_2and8)    ;
% System 3
S3_combined_imp_1_2 = unknownSystem3(imp_1and2)    ;
S3_combined_imp_1_4 = unknownSystem3(imp_1and4)    ;
S3_combined_imp_2_4 = unknownSystem3(imp_2and4)    ;
S3_combined_imp_2_8 = unknownSystem3(imp_2and8)    ;

% Now we can combine the responses AFTER the inputs have been fed through
% the system

% System 1
S1_ind_combined_1and2 = imp_resp_1_1 + imp_resp_1_2     ;
S1_ind_combined_1and4 = imp_resp_1_1 + imp_resp_1_4     ;
S1_ind_combined_2and4 = imp_resp_1_2 + imp_resp_1_4     ;
S1_ind_combined_2and8 = imp_resp_1_2 + imp_resp_1_8     ;
% System 2
S2_ind_combined_1and2 = imp_resp_2_1 + imp_resp_2_2     ;
S2_ind_combined_1and4 = imp_resp_2_1 + imp_resp_2_4     ;
S2_ind_combined_2and4 = imp_resp_2_2 + imp_resp_2_4     ;
S2_ind_combined_2and8 = imp_resp_2_2 + imp_resp_2_8     ;
% System 3
S3_ind_combined_1and2 = imp_resp_3_1 + imp_resp_3_2     ;
S3_ind_combined_1and4 = imp_resp_3_1 + imp_resp_3_4     ;
S3_ind_combined_2and4 = imp_resp_3_2 + imp_resp_3_4     ;
S3_ind_combined_2and8 = imp_resp_3_2 + imp_resp_3_8     ;
% Finally, we will take the responses to combined inputs and see if they
% are equal to the combined responses of their individual outputs

%% System 1 Additivity Check

% We do this by taking the difference between the two 

S1_Additivity_Check_1 = [(S1_ind_combined_1and2 - S1_combined_imp_1_2) < 1e-8]'
S1_Additivity_Check_2 = [(S1_ind_combined_1and4 - S1_combined_imp_1_4) < 1e-8]'
S1_Additivity_Check_3 = [(S1_ind_combined_2and4- S1_combined_imp_2_4) < 1e-8]'
S1_Additivity_Check_4 = [(S1_ind_combined_2and8- S1_combined_imp_2_8) < 1e-8]'

System_1_AdditivityCheck = {'We can see that system 1 is NOT LINEAR because the difference'
    'between the pre-combined and post-combined responses should be zero. Instead, we find here that they.'
    'are different at all points. This tells us that the response can''t be linear as it does not follow the'
    'law of additivity'} 


%% System 1 Additivity Plot Example

% Here I plot one comparison as an eample for good measure. You can see in
% the Y axis that these values differ.

figure
subplot(2,1,1)
stem(S1_combined_imp_1_2, 'r', 'filled')        ;
title('Inputs Combined')                        ;
ylabel('Mean Spike Count')                      ;
xlabel('Time')                                  ;
subplot(2,1,2)
stem(S1_ind_combined_1and2, 'b', 'filled')      ;
title('Individual Responses Combined')          ;
ylabel('Mean Spike Count')                      ;
xlabel('Time')                                  ;
sgtitle('System 1 Linearity Check Responses')   ;

%% System 2 Additivity Check

% Check the difference between the two for system 2
S2_additivity_Check_1 = [(S2_ind_combined_1and2 - S2_combined_imp_1_2) < 1e-8]'
S2_additivity_Check_2 = [(S2_ind_combined_1and4 - S2_combined_imp_1_4) < 1e-8]'
S2_additivity_Check_3 = [(S2_ind_combined_2and4- S2_combined_imp_2_4) < 1e-8]'
S2_additivity_Check_4 = [(S2_ind_combined_2and8- S2_combined_imp_2_8) < 1e-8]'

System_2_AdditivityCheck = {'We can see that system 2 passes our linearity check of additivity.'
    'Differences between the pre-combined and post-combined responses are all zero.'} 


%% System 2 Additivity Plot Example

% I plot an example to see this visually

figure
subplot(2,1,1)
stem(S2_combined_imp_1_2, 'r', 'filled')        ;
title('Inputs Combined')                        ;
ylabel('Mean Spike Count')                      ;
xlabel('Time')                                  ;
subplot(2,1,2)
stem(S2_ind_combined_1and2, 'b', 'filled')      ;
title('Individual Responses Combined')          ;
ylabel('Mean Spike Count')                      ;
xlabel('Time')                                  ;
sgtitle('System 2 Linearity Check Responses')   ;

%% System 3 Linear Check
S3_Additivity_Check_1 = [(S3_ind_combined_1and2 - S3_combined_imp_1_2) < 1e-8]'
S3_Additivity_Check_2 = [(S3_ind_combined_1and4 - S3_combined_imp_1_4) < 1e-8]'
S3_Additivity_Check_3 = [(S3_ind_combined_2and4- S3_combined_imp_2_4) < 1e-8]'
S3_Additivity_Check_4 = [(S3_ind_combined_2and8- S3_combined_imp_2_8) < 1e-8]'

System_3_AdditivityCheck = {'We can see that system 3 also passes our linearity check of additivity'
    'between the pre-combined and post-combined responses are all zero.'} 


%% System 3 Additivity Plot Example

% Lets take a gander at an example comparison
figure
subplot(2,1,1)
stem(S3_combined_imp_1_2, 'r', 'filled')        ;
title('Inputs Combined')                        ;
ylabel('Mean Spike Count')                      ;
xlabel('Time')                                  ;
subplot(2,1,2)
stem(S3_ind_combined_1and2, 'b', 'filled')      ;
title('Individual Responses Combined')          ;
ylabel('Mean Spike Count')                      ;
xlabel('Time')                                  ;
sgtitle('System 3 Linearity Check Responses')   ;

%% pt. 3  Boundary Handling
% Use different n values to determine how the system handles inputs near
% the boundaries. 

% Create impulses towards the "end" of the input vector
impulse_62 = circshift(impulse, 61)  ;
impulse_63 = circshift(impulse, 62)  ;
impulse_64 = circshift(impulse, 63)  ;

%% Create Impulse Responses for System 1 Boundary Handling
imp_resp_1_62 = unknownSystem1(impulse_62) ;
imp_resp_1_63 = unknownSystem1(impulse_63) ;
imp_resp_1_64 = unknownSystem1(impulse_64) ;

%% Plot System 1 Boundary Handling
figure ;
subplot(3,1,1)
stem(imp_resp_1_62, 'k', 'filled')              ;
hold on
stem(62,imp_resp_1_62(62), 'r', 'filled')       ;
stem(63,imp_resp_1_62(63), 'b', 'filled')       ;
stem(64,imp_resp_1_62(64), 'g', 'filled')       ;
title('Position 62 Impulse')                    ;
ylabel('Mean Spike Count')                      ;
xlabel('Time')                                  ;

subplot(3,1,2)
stem(imp_resp_1_63, 'k', 'filled')              ;
hold on
stem(63,imp_resp_1_63(63), 'r', 'filled')       ;
stem(64,imp_resp_1_63(64), 'b', 'filled')       ;
stem(1,imp_resp_1_63(1), 'g', 'filled')         ;
title('Position 63 Impulse')                    ;
ylabel('Mean Spike Count')                      ;
xlabel('Time')                                  ;

subplot(3,1,3)
stem(imp_resp_1_64, 'k', 'filled')              ;
hold on
stem(64,imp_resp_1_64(64), 'r', 'filled')       ;
stem(1,imp_resp_1_64(1), 'b', 'filled')         ;
stem(2,imp_resp_1_64(2), 'g', 'filled')         ;
title('Position 64 Impulse')                    ;
ylabel('Mean Spike Count')                      ;
xlabel('Time')                                  ;
sgtitle('System 1 Boundary Handling: Circular') ;

System_1_BoundaryCheck = {'We can see that system 1 seems to handle boundaries by'
    'in a circular fashion. When an impulse is given at the end of the'
    'vector, we see that the "cut-off" response repeats again at the beginning.'
    'This acts as if these ends are connected, continuing the response at the beginning.'} 


%% Create Impulse Responses for System 2 Boundary Handling
imp_resp_2_62 = unknownSystem2(impulse_62) ;
imp_resp_2_63 = unknownSystem2(impulse_63) ;
imp_resp_2_64 = unknownSystem2(impulse_64) ;

%% Plot System 2 Boundary Handling
figure ;
subplot(3,1,1)
stem(imp_resp_2_62, 'k', 'filled')                  ;
hold on
stem(62:64,imp_resp_2_62(62:64), 'b', 'filled')     ;
stem(1:14,imp_resp_2_62(1:14), 'b', 'filled')       ;
title('Position 62')                                ;
ylabel('Mean Spike Count')                          ;
xlabel('Time')                                      ;

subplot(3,1,2)
stem(imp_resp_2_63, 'k', 'filled')                  ;
hold on
stem(63:64,imp_resp_2_63(63:64), 'b', 'filled')     ;
stem(1:15,imp_resp_2_63(1:15), 'b', 'filled')       ;
title('Position 63')                                ;
ylabel('Mean Spike Count')                          ;
xlabel('Time')                                      ;

subplot(3,1,3)
stem(imp_resp_2_64, 'k', 'filled')                  ;
hold on
stem(64,imp_resp_2_64(64), 'b', 'filled')           ;
stem(1:16,imp_resp_2_64(1:16), 'b', 'filled')       ;
title('Position 64')                                ;
ylabel('Mean Spike Count')                          ;
xlabel('Time')                                      ;
sgtitle('System 2 Boundary Handling: Circular')     ;

System_2_BoundaryCheck = {'System 2 also seems to handle boundaries in a circular fashion.'} 


%% Create Impulse Responses for System 3 Boundary Handling
imp_resp_3_62 = unknownSystem3(impulse_62) ;
imp_resp_3_63 = unknownSystem3(impulse_63) ;
imp_resp_3_64 = unknownSystem3(impulse_64) ;

%% Plot System 3 Boundary Handling

figure ;
subplot(3,1,1)
stem(imp_resp_3_62, 'k', 'filled')          ;
hold on
stem(62,imp_resp_3_62(62), 'r', 'filled')   ;
stem(63,imp_resp_3_62(63), 'b', 'filled')   ;
stem(64,imp_resp_3_62(64), 'g', 'filled')   ;
stem(1,imp_resp_3_62(1), 'o', 'filled')     ;
title('Position 62')                        ;
ylabel('Mean Spike Count')                  ;
xlabel('Time')                              ;

subplot(3,1,2)
stem(imp_resp_3_63, 'k', 'filled')          ;
hold on
stem(63,imp_resp_3_63(63), 'r', 'filled')   ;
stem(64,imp_resp_3_63(64), 'b', 'filled')   ;
stem(1,imp_resp_3_63(1), 'g', 'filled')     ;
stem(2,imp_resp_3_63(2), 'o', 'filled')     ;
title('Position 63')                        ;
ylabel('Mean Spike Count')                  ;
xlabel('Time')                              ;

subplot(3,1,3)
stem(imp_resp_3_64, 'k', 'filled')          ;
hold on
stem(64,imp_resp_3_64(64), 'r', 'filled')   ;
stem(1,imp_resp_3_64(1), 'b', 'filled')     ;
stem(2,imp_resp_3_64(2), 'g', 'filled')     ;
stem(3,imp_resp_3_64(3), 'o', 'filled')     ;
title('Position 64')                        ;
ylabel('Mean Spike Count')                  ;
xlabel('Time')                              ;
sgtitle('System 3 Boundary Handling: Circular');

System_3_BoundaryCheck = {'Finally, system 3 also seems to handle boundaries in a circular fashion.'} 


%% B) Check whether the outputs are sinusoids of the same frequency

% Examine the response of the system to sinusoids with frequencies:
    % - 2(pi)/N         % 4(pi)/N
    % - 8(pi)/N         % 16(pi)/N
    % - AND Random phases

% ** This means that we need to verify the output vector lies completely in
% the subspace containing all the sinusoids of that frequency
% NOTE: Make the input stimuli positive, by adding one to each sinusoid,
% and the responses should then be positive (mean spike counts).

N = 64          ; % # of points
k = [1,2,4,8]   ; % utilized to create the frequency
x = 1:N         ; % domain to plot

% First I create sinusoids for each frequency providing a random phase
% shift
sig1 = cos(2*pi*k(1)/N * x + rand)'     ;
sig2 = cos(2*pi*k(2)/N * x + rand)'     ;
sig3 = cos(2*pi*k(3)/N * x + rand)'     ;
sig4 = cos(2*pi*k(4)/N * x + rand)'     ;

% Then I find their response to the system. I only test System two
% because it was the only system that was found to be an LSI system
sig1_resp = unknownSystem2(sig1) ;
sig2_resp = unknownSystem2(sig2) ;
sig3_resp = unknownSystem2(sig3) ;
sig4_resp = unknownSystem2(sig4) ;

% First we check visually whether or not the outputs are the same frequency
% as the inputs by plotting them...

figure ;
% Signal 1 vs. it's response
subplot(2,4,1)
plot(sig1, 'b')                                     ;
title('Original: 2(pi)/N')                          ;
ylabel('Mean Spike Count')                          ;
xlabel('Time')                                      ;
subplot(2,4,5)
plot(sig1_resp, 'b--')                              ;
title('System Response: 2(pi)/N')                   ;
ylabel('Mean Spike Count')                          ;
xlabel('Time')                                      ;

% Signal 2 vs. it's response
subplot(2,4,2)
plot(sig2, 'r')                                     ;
title('Original: 4(pi)/N')                          ;
ylabel('Mean Spike Count')                          ;
xlabel('Time')                                      ;
subplot(2,4,6)
plot(sig2_resp, 'r--')                              ;
title('System Response: 4(pi)/N')                   ;
ylabel('Mean Spike Count')                          ;
xlabel('Time')                                      ;

% Signal 3 vs. it's response
subplot(2,4,3)
plot(sig3, 'g')                                     ;
title('Original: 8(pi)/N')                          ;
ylabel('Mean Spike Count')                          ;
xlabel('Time')                                      ;
subplot(2,4,7)
plot(sig3_resp, 'g--')                              ;
title('System Response: 8(pi)/N')                   ;
ylabel('Mean Spike Count')                          ;
xlabel('Time')                                      ;

% Signal 4 vs. it's response
subplot(2,4,4)
plot(sig4, 'k')                                     ;
title('Original: 16(pi)/N')                         ;
ylabel('Mean Spike Count')                          ;
xlabel('Time')                                      ;
subplot(2,4,8)
plot(sig4_resp, 'k--')                              ;
title('System Response: 16(pi)/N')                  ;
ylabel('Mean Spike Count')                          ;
xlabel('Time')                                      ;
sgtitle('Original Sinusoids vs. System Output')     ;

Visual_Sinusoid_FrequencyCheck = {'We can get a visual feel for the fact that the frequencies between'
    'the original sinusoids and their responses to system 2 return sinuisoids of the same frequency'
    'and which have been shifted and scaled - as LSI systems should.'} 


%% Checking that the Output Lies Within the Proper Frequency Subspace
% In order to verify this quantitatively (by checking that the output falls
% within the same subspace) we can take the fourier transformation of outputs 
% from System 2 and plot the amplitudes...


% Now take each of their amplitudes
sig1_resp_amp = abs(fft(sig1_resp)) ;
sig2_resp_amp = abs(fft(sig2_resp)) ;
sig3_resp_amp = abs(fft(sig3_resp)) ;
sig4_resp_amp = abs(fft(sig4_resp)) ;

% We only plot the first half of the Fourier Transform because all elements
% in the "second half" are complex conjugates of the first half, so
% it doesn't add any new information.

xplots = linspace(1,31,31) ;

figure ;
subplot(2,2,1)
stem(xplots, sig1_resp_amp(2:(N/2)), 'b', 'filled')
hold on
stem(xplots(1), sig1_resp_amp(2), 'r', 'filled')
title('Signal 1: 2pi/N')
xlabel('Freq.')
ylabel('Amplitude')

subplot(2,2,2)
stem(xplots, sig2_resp_amp(2:(N/2)), 'b', 'filled')
hold on
stem(xplots(2), sig2_resp_amp(3), 'r', 'filled')
title('Signal 2: 4pi/N')
xlabel('Freq.')
ylabel('Amplitude')

subplot(2,2,3)
stem(xplots, sig3_resp_amp(2:(N/2)), 'b', 'filled')
hold on
stem(xplots(4), sig3_resp_amp(5), 'r', 'filled')
title('Signal 3: 8pi/N')
xlabel('Freq.')
ylabel('Amplitude')

subplot(2,2,4)
stem(xplots, sig4_resp_amp(2:(N/2)), 'b', 'filled')
hold on
stem(xplots(8), sig4_resp_amp(9), 'r', 'filled')
title('Signal 4: 16pi/N')
xlabel('Freq.')
ylabel('Amplitude')
sgtitle('Amplitudes for Each Sinusoid')

Quant_Sinusoid_FrequencyCheck = {'We can see in these plots that the only position which elicits an amplitude response'
    'is the position which matches the frequency. Since we know that each component of the Fourier Transformation'
    'represents the amount of amplitude response that a system returns given a certain frequency, we can see that these match.'
    'Responses for sinusoidal frequency: '
    '2(pi)/N --> frequency of 1 (b/c we factor out the 2*pi)'
    '4(pi)/N --> frequency of 2'
    '8(pi)/N --> frequency of 4'
    '16(pi)/N --> frequency of 8'
    'We can see that each sinusoid"s amplitude response takes place in the cooresponding position of'
    'the fourier transformed response. This confirms that the Output Lies Within the Proper Frequency Subspace'} 

%% C) Lets talk about Linear Shift Invariance

%% The Change in Amp. and Phase from Input to Output in FT Space

% For this portion, I'm going to utilize only sig1 from
% earlier

% First I find the fourier transformed amplitudes and phases of the input
signal1_input_amps = abs(fft(sig1)) ;
signal1_input_phases = angle(fft(sig1)) ;

% Now take the phase - amplitude of sig1_resp already found above
sig1_resp_phase = angle(fft(sig1_resp)) ;

% Take the index, and then drop the duplicate
sig1_max_amp_index = find(sig1_resp_amp == max(sig1_resp_amp)) ;
sig1_max_amp_index(2) = [] ;

% We find the ratio of change between the input and output at the most
% responsive frequency
sig1_max_amp_differences_ratio = sig1_resp_amp(sig1_max_amp_index)/ signal1_input_amps(sig1_max_amp_index) 
sig1_phase_differences = sig1_resp_phase(sig1_max_amp_index) - signal1_input_phases(sig1_max_amp_index)  ;

% Now we create an impulse response, put it through unknownSystem2, and
% then see that these differences match the changes in the above phase and
% amplitudes

impulse = zeros(64,1)   ;
impulse(1) = 1          ;
impulse_response = unknownSystem2(impulse) ;

% Take the fft amplitude of the original impulse and it's response
impusle_input_amps = abs(fft(impulse)) ;
impusle_response_amps = abs(fft(impulse_response)) ;

% Take the fft phases of the original impulse and it's response
impulse_input_phases = angle(fft(impulse)) ;
impulse_output_phases = angle(fft(impulse_response)) ;

impulse_amp_difference_ratio = impusle_response_amps(sig1_max_amp_index) / impusle_input_amps(sig1_max_amp_index) ;
impulse_phase_differences = impulse_output_phases(sig1_max_amp_index) - impulse_input_phases(sig1_max_amp_index) ;

% Round a few spots back so we can get a clean comparison
impulse_amp_difference_ratio = round(impulse_amp_difference_ratio,4)          ;
impulse_phase_differences = round(impulse_phase_differences,4)
sig1_max_amp_difference_ratio = round(sig1_max_amp_differences_ratio,4)
sig1_phase_differences = round(sig1_phase_differences,4)

% Test that these are equal
amp_diff_check = isequal(sig1_max_amp_difference_ratio,impulse_amp_difference_ratio)
phase_diff_check = isequal(impulse_phase_differences,sig1_phase_differences)

Answer = {'Here we can see that the change in phase and amplitude within the fourier'
            'frequency domain between the impulse vector and it''s response predicts'
            'the changes that we also saw above with the sinusoid sig1.'}

%% Testing how Scaling the Sinusoid Changes Amplitude in the Fourier Frequency Domain
% We'll only be checking this for unknownSystem2 because it is the only
% system which has passed the previous tests. Here we will scale and shift
% an input signal and see if we can see those changes in the fourier
% frequency domain. I did this before doing anything else and I couldn't
% bring myself to delete it.

N = 64              ;
scale_factor = 2    ;
shift = 10          ;
test_x = -N/2:N/2-1 ;

test_signal = cos(2*pi*1/64 * test_x )'             ;
test_shifted = cos(2*pi*1/64 * test_x + shift )'    ;
test_scaled_and_shifted = test_shifted*scale_factor ;

test_sig_response = unknownSystem2(test_signal)                             ;
test_sig_shiftANDscale_response = unknownSystem2(test_scaled_and_shifted)   ;

test_sig_norm_amps = abs(fftshift(fft(test_sig_response)))      ;
test_sig_norm_phase = angle(fftshift(fft(test_sig_response)))   ;

test_sig_shiftANDscaled_amps = abs(fftshift(fft(test_sig_shiftANDscale_response)))      ;
test_sig_shiftANDscaled_phase = angle(fftshift(fft(test_sig_shiftANDscale_response)))    ;

figure ;
subplot(2,1,1)
stem(test_x, test_sig_norm_amps, 'b', 'filled')
title('Original Test Signal: Amplitudes')
xlabel('Freq.')
ylabel('Amplitude')

subplot(2,1,2)
stem(test_x , test_sig_shiftANDscaled_amps, 'b', 'filled')
title('Scaled Test Signal (by 2): Amplitudes')
xlabel('Freq.')
ylabel('Amplitude')

% We can see visually that the maximum amplitude factor doubles with the 

% We can also check the amplitudes scaling factor
max_test_scaled_amp = max(test_sig_shiftANDscaled_amps) ;
max_test_norm_amp = max(test_sig_norm_amps)             ;
scale_test =  max_test_norm_amp * scale_factor == max_test_scaled_amp  

amp_difference = max(test_sig_norm_amps) - max(test_sig_shiftANDscaled_amps)

%% Linear Shift Invariant Guarantee?
LSI_Guarantee_Answer = {'No, this doesn"t GUARANTEE either linearity or shift invariance'
    'The simplest reason is that all "useful" linear/shift invariant systems are technically'
    'bounded. They will all break if you expose them to an extreme situation. For example, '
    'A seemingly linear system will certainly break if you put it "close to the sun". Additionally,'
    'We can never know that a system is shift invariant because we will never be able to test'
    'all "possible shifts". For example, if we use time as a domain to shift along, we will never'
    'be able to predict how something will respond in the future - at least without a time machine ;).'
    'For the other systems, I address which tests they fail earlier in this script as I test them.'}


