% In general, the whole homework relies on one principle idea. When you
% take the fft() of a signal, and plot the amplitudes. What you see is the
% amount of amplitude that will be returned at each frequency/period.
% (These are inverses of one another). So if you see a response at x = 3,
% then that is the sampling rate (3 cycles/unit measured) that gives you a
% response.

% Question 1

% (C) For this one you want to create an impulse response that has a value
% of 1 at the frequency location that you found in part B to elicit a
% resposne. Then you translate this into the fft() space and plot the
% amplitudes and phases. I actually haven't plotted the phases yet so I'm
% not sure what they will look like.

% Question 2

% I found the actual question for this to be wildly confusing and to me
% seems to be CLEARLY asking us to do something different from what Ionotan
% told me we need to do.

% (A) for this, you create eight 8-D impulse vectors (so all zeros except
% for 1) in each position. Then you convolve each with r and the put them
% together into one matrix. You are creating the diagonal looking matrix
% that can be gound on slide 6 of the LSIsystems-slides file ("Convolution
% Matrix") is the title of that slide.

% (B) For this one you can look at the documentation for the conv()
% function. It describes the different ways in which the function can
% handle the boundaries. This has to do with whether or not it will zero
% pad or not.

% (C) For this you just create a single cycle 48 vector. To do this you
% can do cos(2*pi/N*x) where N = 48. (Can plot to check this works). If you
% plot the convolved response you will see that the edges are weird.
% Basically the default boundary handling using zero padding, in order to
% fix this, you can pass the 'valid' argument into the conv() function and
% then replot.

% Everything below here I wrote up super quick and may not make total
% sense, let me know if you have any questions

% question 3

% Period = inverse of frequency - cycles per sample

% A - Create the gabor filter
% use fftshift to center it
% This is a gabor filter
% Google gabor filters and talk about why they look the way they look

% B - To find the largest response, you need to take the DFT of the gabor
% filter and then plot the amplitudes, find the max value, and then take
% that frequency. To find the 25% frequncy, you find the frequency who gets
% a amplitude value that is 25% of the max amplitude frequency

% C - Create the frequencies that you found in the previous part of the
% problem and then convolve them, take the DFT and then confirm that the
% amplitude response is what you'd expect. They wont be exaclty right b/c
% we are dealing with discrete sinusoids and not continuous ones.

% Question 4:

% A - Create the your own convolution function. You should be creating the
% same matrix that you created in the second question. Confirm that it
% works by providng vectors of the length mentioned in the problem.
% Visualize as an image. it should look like the visual diagnol form of the
% actual matrix

% B - Do a least squares regression by taking the SVD of X and then working
% into the star star space. Then plot it

% C - Take the fourier transform of the HRF, then convert to amplitude,
% then plot them after fftshifting to center on zero. Describe this type of
% filter. 

% Question 5: 

% A - Take every fourth value of the sig vector. Plot both of these vectors on 
% top of each other. Talk about how this works as a subsampling.
% - linear, because you are sampling at equal distances
% - not shift invariant because you if you start sampling at 4 vs. 5 you
% will get a different end results

% B - Take the Fourier Transformation of sig and then plot the amplitudes
% after fftshifting them around the center. Find which amplitude frequency
% returns the highest value and then CONVERT THIS TO Hz BY MULTIPLYING BY
% 100Hz

% C - Sample at every second, third, and fourth point - fill the other
% values with zeros. What is the relationship for these plots. Why do we
% see that certain peak values are higher than the previous original
% amplitude plot's peak values. This is where you need to look up nyquist
% samppling and "sample folding"















