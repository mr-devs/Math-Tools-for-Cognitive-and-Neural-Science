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

% What will make it work?

%% Projection practice

% define two vectors
aa = [10; 10];
bb = [20; 9];

%what is their dot product? Do NOT use dot()

%project aa onto the line defined by bb
%what is its length? what is its direction?

%% 2D --> 1D projection
% You are tasked with finding a new 'purple' for NYU! The current purple
% has [79, 135] (Math on board) 
load('test_purple.mat')
nyu_purple = [79, 135]'; %alternatively, [79;135]
unit_purple= []; %define the purple axis!

% how purple is NYU's purple?
mag_nyu = []; %magnitude of NYU vector's projection onto purple axis;
% how 

ortho_purple = []; % How do we find the orthogonal vector? 

x_purple = [97; 111]; 
mag_x = []; %how purple is this x_purple?

% Load test_purple.mat into Workspace
% This matrix contains a 2x10 matrix comprising ten (R,B) column vectors. 


%FORLOOP START


%FORLOOP END



%can we do this all in one shot and avoid loops? (hint: matrix multiply)


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
% Recall how we used the 

purple_transform = unit_purple';

u = 1; %1x1
s = [1, 0]; %1x2
vt = [unit_purple, ortho_purple]'; % 2x2

disp('purple_transform is unit_purple transposed')
% purple_transform

disp('what is this chain of linear transformations?')
% u*s*vt
