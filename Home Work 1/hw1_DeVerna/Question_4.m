%% 4 A Simple Visual Neuron


%% Question
%{
Suppose a retinal neuron in a particular species of toad generates
responses that are a weighted sum of the (positive-valued) intensities of light
that is sensed at 6 localized regions of the retina. The weight vector is 
[1, 3, 8, 8, 3, 1]. 
%}

weighting_vector = [1;3;8;8;3;1] ;

%% (a) Is this system linear? If so, how do you know? If not, provide a counter-example.
%{ 

This system must be linear because the system's operations (a weighted sum) are
multiplication and addition. To be linear, a system must obey the rules
of superposition - homogeneity and additivity - both multiplication and
addition are allowed within a linear system.

We can also prove this by creating two different input vectors and testing
whether or not they obey the laws of super position.

By creating a second input which is input_stumli1 scaled by 2, we can see
if the response obeys the laws of homoegeneity.

%}

input_stimuli1 = [1,0,0,0,0,0] ;
input_stimuli2 = [2,0,0,0,0,0] ;


% Because input_stimuli2 is input_stimuli1 scaled by 2, putting these
% inputs into the system should return outputs which are also scaled by two
% - if the system is linear.

output1 = sum((input_stimuli1 * 2) * weighting_vector)
output2 = sum((input_stimuli2 * 2) * weighting_vector)

% Now we can check if output1 * 2 = output 2

if output1 * 2 == output2
    disp('The system is linear!')
else
    disp('The system is NOT linear!')
end

%% (b) What unit-length stimulus vector elicits the largest response in the neuron?
% Explain how you arrived at your answer.

%{
If we look at this problem as getting the magnitude of a projection onto a
unit vector we can see quicky that the larger the ANGLE between the
weighting vector and the unit vector, the smaller a response we will get
from the frog's lovely neuron. Thus, we should take input the unit_vector
that is ON TOP OF the weighting vector itself. The vector going in the
other direction would not physically realizable (it would be negative).

Said another way, we should simply find the weighting vectors unit vector
by dividing by it's length.
%}

length_weighting_vector = sqrt(weighting_vector'*weighting_vector) ;

unit_weight_vector = weighting_vector  / length_weighting_vector ;

largestNeuronResponse = weighting_vector'*unit_weight_vec ;

answer = 'The largest response that we can get from this frog''s neuron = %f.' ;

sprintf(answer,largestNeuronResponse)

%% (c) What physically-realizable unit-length stimulus vector produces the smallest response in this neuron?
% - Explain. 
%{
Since this linear system takes the following form:

response = w1*s1 + w2*s2 + ... + w_n*s_n

where:
w = the weighting vector and
s = the stimulus input

this takes the more concrete form of:

response = 1*s1 + 3*s2 + 8*s_3 + 8*s_4 + 3*s_5 + 1*s_6

As alluded to in part (b) taking a unit vector which creates the largest
angle between the weighting vector and the unit vector gives you the
SMALLEST response. This is similar to projecting onto an orthogonal unit
vector.

Given the form of the response system, we can cycle through orthogonal
bases and clearly see that input_stumulus1 and input_stimulus6 gives us
identical answers, but they are the lowest.

%}

input_stimuli1 = [1,0,0,0,0,0] ;
input_stimuli2 = [0,1,0,0,0,0] ;
input_stimuli3 = [0,0,1,0,0,0] ;
input_stimuli4 = [0,0,0,1,0,0] ;
input_stimuli5 = [0,0,0,0,1,0] ;
input_stimuli6 = [0,0,0,0,0,1] ;

output1 = input_stimuli1 * weighting_vector
output2 = input_stimuli2 * weighting_vector
output3 = input_stimuli3 * weighting_vector
output4 = input_stimuli4 * weighting_vector
output5 = input_stimuli5 * weighting_vector
output6 = input_stimuli6 * weighting_vector

answer = 'One of the two input vectors which gives the lowest answer is input_stimulus6'

sprintf(answer)