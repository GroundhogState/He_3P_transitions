power_data = [0.10 	3.31	744396453.80 .25 	0.27    .05     4.05  .84   12508
    0.15 	5.14 	744396454.18	.170	0.53    .07     5.08  0.71  19248
    0.20 	6.5 	744396454.67	.13 	0.42    .05     3.57  .42   18604
    0.25 	8.9 	744396453.03	.17 	0.64    .07     5.31  .57   23202];
%   V       mW      freq            unc     ratio  rat unc  width w unc num 
% Recorded PD ranges: [0.0475,0.06450.0817,0.099] resp, with '0' reading 0.011
pd_range = [0.0475,0.0645,0.0817,0.099]-0.011;
wm_offsets = [-0.030973,1.519872]; %F3 and F4 respectively
waist = 4.1; %cm
voltage = power_data(:,1);
% power=power_data(:,2);
power = [3.316,4.933,6.496, 8.047]';
power_err = [0.059;.043;.044;.039];
freq_diff=power_data(:,3)-mean(power_data(:,3));
f_unc=power_data(:,4);
ratio = power_data(:,5);
ratio_err = power_data(:,6);
number = power_data(:,7);
intensity = 2*power./(pi*waist^2);

nfit = 4;
power_slope=mean(freq_diff(1:nfit)./power(1:nfit));
p_slope_err = std(freq_diff(1:nfit)./power(1:nfit))/2;
ratio_X = [ones(size(power)),power];
ratio_mdl = ratio_X\ratio;

plot_X = linspace(0,10,2);
ratio_slope=mean(ratio(1:nfit)./power(1:nfit));
slope_err = std(ratio(1:nfit)./power(1:nfit))/2;

f1=stfig('Intensity shift');
clf;
% subplot(2,1,1)
% hold on
% errorbar(power,freq_diff-mean(freq_diff),f_unc,'kx')
% plot([0,10],[0,10*power_slope],'k')
% plot([0,10],freq_mdl(2)*[0,10]+freq_mdl(1))
% xlim([0,10])
% box off
% % xlabel('Beam power (mW)')
% ylabel(sprintf('Centre freq. - %.1f (MHz)',mean(power_data(:,3))))
% set(gca,'FontSize',12)
% 
% subplot(2,1,2)
hold on
plot([0;power(1:nfit)], [0;ratio_slope*power(1:nfit)],'k')
plot([0,10],[[1,0];[1,10]]*ratio_mdl)
errorbar(power(1:nfit),ratio(1:nfit),ratio_err,ratio_err,power_err,power_err,'k.','MarkerSize',10)
% plot([0,10],ratio_mdl(2)*[0,10]+ratio_mdl(1))
xlim([0,10])
xlabel('Beam power (mW)')
ylabel('Peak N$_{\textrm{loss}}$/N$_{\textrm{total}}$')
set(gca,'FontSize',13)
% title('$2^3P_2 - 5^3D_1$ Signal linearity')
legend('Linear fit','Affine fit','Data','Location','SouthEast')
suptitle('$2^3P_2 - 5^3D_1$ signal versus beam power')
% 
% 
% stfig('loose')
% plot(freq-mean(freq),ratio,'kx')
% xlabel('Freq - mean')
% ylabel('Ratio signal')
