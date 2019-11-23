%% HW # 4 - Question 1 - Math Tools - Matthew DeVerna

clear all
close all

%% Question 1: Bayes Rule and Eye Color
% A male and female chimpanzee have blue and brown eyes, respectively. The
% brown-eyed allele can be denoted as a capital B, where as the blue-eyes
% allele can be represented as a lowercase b. Assume a simple genertic
% model in which the gene for brown eyes is always dominant (so that the
% trait of blue eyes can only arise from two blue eyed genes, but the trait
% of brown eyes can arise from two brown eyed genes, or, one brown and one
% blue gene. 

% BLUE gene = b
% BROWN gene = B

% You can assume:
% 1) the probability of the mother being BB is 50% and the probability of
% her being Bb is 50%
% 2) the *a priori* probability that each of the four gene configurations
% is equally possible. 

% For each question, provide the math, and explain your reasoning.

%% (A) Imagine the chimps have a single child with brown eyes.
% Given this situation, what is the probability that the female chimp has a
% blue-eyed gene?

% We know that the mother has the below probabilities for Bb and BB:
mother_Bb = .5 ;
mother_BB = .5 ;

% After creating punnett squares we also know that:
p_brown_child_mom_Bb = .5   ;
p_brown_child_mom_BB = 1    ;

% The questions asks us to find:
% P(mother_Bb | child_brown_eyes)

% This inverts with Bayes Rule to:
% p(child_brown_eyes | mother_Bb) * p(mother_Bb) / p(child_brown_eyes)

% The problem becomes much simpler if you decompose p(child_brown_eyes) into:
% P(p_Child_brown | mother_Bb) * P(mother_Bb) + P(p_Child_brown | mother_BB) * P(mother_BB)

% So really, we already have everything we need except p(child_brown_eyes), so
% we populate the decomposed form of p(child_brown_eyes) and then populate from
% the punnett square probabilities...

numerator = p_brown_child_mom_Bb * mother_Bb                                ;
denom = p_brown_child_mom_Bb * mother_Bb + p_brown_child_mom_BB * mother_BB ;

p_mother_Bb_given_child_brown_eyes = numerator/denom

%% (B) A second child with brown eyes.
% We can now look at this as, what is the probability of getting two
% brown-eyed children given that the mother has the gene configuration Bb?

% In order to address this, we will need to update our priors from the
% p(brown eyed child | mother_Bb) to p(TWO brown eyed children | mother
% Bb).

% We find this by multiplying the previous probability by itself, so the
% formula only changes slightly to the below...
numerator = p_brown_child_mom_Bb^2 * mother_Bb                              ;
denom = (p_brown_child_mom_Bb * mother_Bb + p_brown_child_mom_BB * mother_BB)^2 ;

p_mother_Bb_given_TWO_child_brown_eyes = numerator/denom

%% (C) Creating a probability function 
% In order to generalize, suppose the chimps have N children with brown
% eyes - express the probability as a function of N.

% To generalize this, we realize that we multiply p_brown_child_mom_Bb by
% itself N times. Thus, we can just raise this probability to N and it will
% work for any N.

% Here I plug in 1000 for N, which is simulating the idea that the mother
% would have 1000 children (yikes...) all with brown eyes. The probability
% of her having a blue-eyed gene, while virtually never having a blue eyed
% child, EVEN THOUGH the father has only blue-eyed genes is extremely rare.

% # of children this baby factor chimp is having, and all in the name of probability...
N = 1000; 

numerator = p_brown_child_mom_Bb^N * mother_Bb                              ;
denom = (p_brown_child_mom_Bb * mother_Bb + p_brown_child_mom_BB * mother_BB)^N ;

p_mother_Bb_given_N_child_brown_eyes = numerator/denom

