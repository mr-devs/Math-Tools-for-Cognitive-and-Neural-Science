%% Lab 3
% This lab will teach you and give insight on projections, linear
% transformations.


%% Matrix multiplication
% Multiply these two matrices. 
% Remember (m_rows x n_cols)(n_rows x p_cols) = (m_rows x p_cols)
 
A = [1 4 3 2;
    3 2 8 1];

B = [1 8 1 0;
    6 9 2 3;
    4 2 0 1];

%A*B = []; % Will this work?
% Answer: no; always check sizes of matrices before you multiply, it will
% save you hours of headache down the line. 

% What will make it work?
% Hint: transposes... 

A*B'

B*A'

% Note how the resultant products of these matrix multiplies are transposes
% of each other! 

% (XY)' = Y'X'


%% Projection practice

% define two vectors
aa = [10; 10];
bb = [20; 9];

%what is their dot product? Do NOT use dot()
dotted = aa'*bb; 

len_a = sqrt(aa'*aa);
len_b = sqrt(bb'*bb);

% dotted = len_a * len_b * cos(angle_ab)
% therefore angle_ab = ?
angle_ab = acosd(dotted/(len_a*len_b)) ;

%project aa onto the line defined by bb
%what is its length? what is its direction?

u_b = bb/len_b; %define the unit vector that lies on the axis of proj.
len_proj_ab = aa'*u_b; % length of the projection 'shadow'
proj_ab = len_proj_ab * u_b; % length * direction = vector!

%% 2D --> 1D projection
% You are tasked with finding a new 'purple' for NYU! The current purple
% has [79, 135] (Math on board) 
load('test_purple.mat')
nyu_purple = [79, 135]'; %alternatively, [79;135]
unit_purple= [1; 1]/sqrt(2); %define the purple axis!

% how purple is NYU's purple?
mag_nyu = nyu_purple'*unit_purple; %magnitude of NYU vector's projection onto purple axis;

% How do we find the orthogonal vector? 

ortho_purple = nyu_purple - mag_nyu*unit_purple; %subtract its projection onto the purple axis!
ortho_purple = ortho_purple/sqrt(ortho_purple'*ortho_purple); %normalize

x_purple = [97; 111]; 
mag_x = x_purple'*unit_purple; %how purple is this x_purple?

% Load test_purple.mat into Workspace
% This matrix contains a 2x10 matrix comprising ten (R,B) column vectors. 


%FORLOOP START
proj_test = zeros(10,1);
for ii = 1:10
    proj_test(ii) = test_purple(:,ii)'*unit_purple;
end
proj_test
%FORLOOP END

%can we do this all in one shot and avoid loops? (hint: matrix multiply)
unit_purple'*test_purple %one line!


%% plot the purples
% Why are most of them 150 units of purple despite them being very
% different values?

figure
hold on
plot([0,255],[0,255],'--m','linewidth',2);
quiver(0,0,purple_nyu(1),purple_nyu(2),'k','linewidth',1); %black one is original nyu purple
for ii = 1:10
    
    if ii<10
        col = [0 0 1];
    else
        col = [1 0 0];
    end
    quiver(0,0,test_purple(1,ii),test_purple(2,ii),'color',col)
    
end
axis([0,150,0,150])
axis square
xlabel('Redness')
ylabel('Blueness')

%% 3D --> 2D projection

% data vectors
% [height (in) , weight (lbs), Shoe size (Eur)]
lyndon = [68, 140, 9.5];
person1 = [68, 140, 7]; 
person2 = [68, 140, 6]; 
person3 = [68, 140, 6.5]; 
person4 = [68, 140, 14]; 

% put them all into a matrix
people = [lyndon; person1; person2; person3; person4];

[n_people, n_variables] = size(people);

% quiver3(starts(:,1), starts(:,2), starts(:,3), ends(:,1), ends(:,2), ends(:,3))
close all
figure 
hold on
for pp = 1:n_people
    quiver3(0, 0, 0, people(pp,1), people(pp,2), people(pp,3))
%     plot3([0 people(pp,1)],[0 people(pp,2)],[0 people(pp,3)])
end
xlabel('Height (in)')
ylabel('Weight (lbs)')
zlabel('Shoe size')
view(3)
axis square
grid on
legend('Lyndon','person1','person2','person3','person4')

%% Matrix = linear transformation
% Back to purple space... 
% We can decompose the unit purple vector into a chain of matrix multiplies
% Remember, a scalar and vector are both special cases of matrices.
% Just like how you can write 5 = 1/3 * 5 * 3 , you can decompose linear
% transform matrices into a sequence of matrix multiplies that can help you 
% better understand what the transformation is doing... 

purple_transform = unit_purple';

s = [1, 0]; %1x2
vt = [unit_purple'; ortho_purple']; % 2x2

disp('purple_transform is unit_purple transposed')
% purple_transform

disp('what is this chain of linear transformations?')
% s*vt'
s*vt 

% We will go into more depth about what is happening here during the SVD
% lecture. But for now, note that the matrix s = [1,0] is multiplying the
% matrix vt = [.7071, .7071; -7071 .7071]. This is a (1x2) matrix
% multiplying a (2x2) matrix, so we expect their matrix product to be of
% size (1x2). The first row of vt is unit_purple, and the second row of vt
% is ortho_purple (the unit vector that is orthogonal to unit_purple). What
% s*vt is effectively doing is 'selecting' the first row of vt, and
% ignoring the 2nd row of vt, since s=[1,0]. The product is thus [0.7071,
% 0.7071] which is unit_purple'!

unit_purple' == s*vt

% What we describe here is a very powerful statement, but looks as trivial
% as 5 = 1/3 * 5 * 3. We decomposed unit_purple' into a selector matrix s
% that 'selected' the unit_purple component of vt, and zeroed the
% ortho_purple copmonent of vt. Thus, we project a vector X onto
% unit_purple by unit_purple'*X, we are saying this equivalent to saying
% s*vt*X. In other words, a projection of X onto unit purple, and a vector
% orthogonal to unit_purple, and then zeroing the orthogonal component
% while keeping the component that lies on unit_purple. This concept is the
% heart of the singular value decomposition.
