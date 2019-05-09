data = [0.10 	3.31	744396515.790 .264 	0.36	1.67
    0.15 	5.14 	744396516.282	.180	0.62	2.20
    0.20 	6.5 	744396516.508	.128 	0.74	2.28
    0.25 	8.9 	744396515.106	.157	0.84	2.37];
wm_offsets = [-0.030973,1.519872]; %F3 and F4 respectively
waist = 4.1; %cm
p=data(:,2);
f=data(:,3);
e=data(:,4);
intensity = 2*p./(pi*waist^2)

sfigure(15852);
clf;
set(gcf,'color','w');
errorbar(p,f-mean(f),e,'kx')
xlabel('Beam power (mW)')
ylabel(sprintf('Centre freq. - %.2f (MHz)',mean(f)))
title('2^3P_2\rightarrow 5^3D_1 gap versus beam power')
lande_sg