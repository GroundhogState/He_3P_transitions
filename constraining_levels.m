% Lost the doppler width script :(
disp('Calculating:')

h=6.62607015e-34;
c=299792458;
% Well, Lach and Pachucki (2001) give the M1 transition at 
delta_E = 0.1065403108; %atomic units!!!
A = 6.484690e-9; %Hz
hartree = 4.3597447222071e-18; %J, atomic unit conversion
E_lach = delta_E * hartree;
f_lach = E_lach / h
% Derevianko prediction?

% From Drake;
% Consider transition from another state; whatever can get there!
E_upp = 183236.791701; % cm^-1
E_upp_err = 0.000067;% cm^-1
E_low_err = 0.0011; % cm^-1
E_tot_err = E_low_err + E_upp_err;
E_meta=159855.9726;% cm^-1
% metastable transition
l_vac = 625.563e-10; %m, lambda
% So the difference is
K_fbd = 100*(E_upp-E_meta); %m^-1
E_fbd = h*c*K_fbd; %J
format longg
F_fbd = c*K_fbd

diff = (f_lach- F_fbd)/1e9
f_err = c*100*E_tot_err/1e9
%% Drake's big table
E_23S1 = 159855.974330;
E_23P2 = 169086.766473;

%% Working out constraints we provided
f_meas = [727303247,4;
	744396515,20;
    744396235,20;
	744396204,20;
	744430345,20];

E_33D1 = [366018892.97,0.02];%Mhz %Ionization energy from Pachucki
E_31D2 = [365917749.02,0.02];%Mhz
Luo_23P0_33D1 = [510059755.352,0.28]; %Mhz
Pachucki_23P0_23P2 = [27599.0572+4309.0742,0.0016];% MHz
% So we want to fix the IE of the measured transitions?
IE = E_33D1 + Luo_23P0_33D1 +  Pachucki_23P0_23P2 - f_meas %MHz






% 'Through connections with the predicted 3^3D and 3^1D levels...'
% So we ASSUME the 3D levels and then constrain the  upper levels
% through calc E_target = E_3D - E_23P2 + h*f_transition
% Or really, E_target = E_23P2 + h*f_transition, with uncertainty from the
% 23P2

