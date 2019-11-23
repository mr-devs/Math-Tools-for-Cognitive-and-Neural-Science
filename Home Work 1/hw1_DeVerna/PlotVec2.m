function PlotVec2(MATrix)
%Function which takes matrix of height two, plotting each column on a 2D
%axis

%{   
MATrix needs to be a matrix with a height of 2. Should this matrix not have
height two, an error will be thrown indicating as such.
%}

if size(MATrix,1) == 2
    for ii = [1: length(MATrix)]
    plot([0,MATrix(1,ii)], [0,MATrix(2,ii)], 'LineWidth',1.5)
    hold on
    plot([MATrix(1,ii)], [MATrix(2,ii)], 'k*' ,'LineWidth',1.5)
    end

title('Plotting of 2D Vectors')
xlabel('x-axis')
ylabel('y-axis')
xlim([-1 1])
ylim([-1 1])
grid on

axis equal

else
    error('Matrix needs to have two rows. Please check the dimensions of your matrix. :)')


end
