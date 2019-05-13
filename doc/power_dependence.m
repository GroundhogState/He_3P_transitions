data = [0.10 	3.31	744396515.790 .264 	0.36	1.67  12508
    0.15 	5.14 	744396516.282	.180	0.62	2.20  19248
    0.20 	6.5 	744396516.508	.128 	0.74	2.28  18604
    0.25 	8.9 	744396515.106	.157	0.84	2.37  23202];
%   V       mW      freq            unc     ratio   width  num 
wm_offsets = [-0.030973,1.519872]; %F3 and F4 respectively
waist = 4.1; %cm
p=data(:,2);
f=data(:,3);
e=data(:,4);
r = data(:,5);
n = data(:,7);
intensity = 2*p./(pi*waist^2);


X = linspace(0,10,100);

f1=sfigure(15852);
clf;
set(gcf,'color','w');
errorbar(p,f-mean(f),e,'kx')
xlabel('Beam power (mW)')
ylabel(sprintf('Centre freq. - %.2f (MHz)',mean(f)))
title('2^3P_2\rightarrow 5^3D_1 gap versus beam power')

f2=sfigure(24582);
clf;
set(gcf,'color','w');
subplot(2,1,1)
plot(p,r,'ko')
hold on
% plot(X,.2*X.^(2/3))
% plot(X,.1*X)
% plot(X,.1*X.^(2/3)+0.05*X)
title('Ratio')
xlim([0,10])
ylim([0,1])
subplot(2,1,2)
plot(p,n,'ko')
hold on
% plot(X,5.5e3*X.^(2/3))
% plot(X,3e3*X)
title('Number')
xlim([0,10])
ylim([0,2.5e4])
xlabel('Power')
suptitle('2^3P_2\rightarrow 5^3D_1 Signal linearity')

[n./p,r./p]
