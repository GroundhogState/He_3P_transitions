%Analysis of 5^3D states
g_state = '2_3P_2_2';
% ITC 1         ITC 2
data ={744396158380000,744396180286000,'5_3D_3_1'
           744396191764000,744396205362000,'5_3D_2_1'
           744396220906000,744396225174000,'5_3D_2_2'
           744396235218000,744396225174000,'5_3D_3_3'
           744396448622000,744396473900000,'5_3D_1_1'};
B_exp = [16.5,10.5];
sfigure(404);
clf
B_plot = B_exp.*ones(size(data,1),2);
for ii =1:size(data,1)
    hold on
    plot(B_plot(ii,:),cell2mat(data(ii,1:2)))
end
B=linspace(0,20);       
e_state = '5_3D_3_1';
lines = zeeman_state2state(B,g_state,e_state,const);

xlabel('Magnetic field strength (G)')
ylabel('Frequency (Hz)')