

%% Complex numbers

disp('Matlab defaults the imaginary number to i and j.')
disp('i=')
disp(i)
disp('j=')
disp(j)
disp('sqrt(-1)=')
disp(sqrt(-1))

%% Real and imaginary parts of complex numbers
cc = 5 + 8j;
disp('c =')
disp(cc)
sprintf('real(c) = %f', real(c))
sprintf('imag(c) = %f', imag(c))

%%
fs = 1000; 
dt = 1/fs; 
stop_time = .1;
xx = (0:dt:stop_time)';

fc =25;

omega = 2*pi*fc;
z1 = cos(omega*x) + 1j*sin(omega*x);
hold on
plot(real(z1))
% plot(imag(z1))



%%
close all
hold on

lims = [0,11,0,2];

subplot(211)
xx=[1,0,0,0,0,0,0,0,0,0];
yy=[0,.1,.25,.5,.75,.5,.3,.2,.1,0];
stem(xx,'filled')
axis(lims)

subplot(212)
stem(yy,'filled','r')
axis(lims)


figure

subplot(211)
x1 = circshift([1,0,0,0,0,0,0,0,0,0],1);
y1 = circshift([0,.1,.25,.5,.75,.5,.3,.2,.1,0],1);
stem(x1,'filled')
axis(lims)

subplot(212)
stem(y1,'filled','r')
axis(lims)

figure
subplot(211)
stem(xx+x1,'filled')
axis(lims)

subplot(212)
stem(yy+y1,'filled','r')
axis(lims)

%%


xx = [1,1,1]';

yy = ones(5,1);

conv(yy,xx,'valid')



