f_exp = [744396159.014 744396190.235 744396221.123 % stage 1
    744396180.504 744396204.973 744396222.922]; % stage 2
f_add = [744396159.014 744396179 744396190.235 744396221.123 % stage 1
    744396180.504 744396190 744396204.973 744396222.922]; % stage 2
f_single = [744396190,;744396204];
w_single = [4.90,4.83];

f_err = [0.21,0.16,0.16;0.14,0.13,0.16];
w_exp = [3.44,8.56,3.67;4.34,7.33,4.75];
B_exp = [18.25,11.43]; %Gauss;

% sfigure(23408)
% clf;
% plot(B,f,'kx')


%% init
% clear all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
addpath(genpath(core_folder));
zeeman_folder = fullfile(fileparts(this_folder),'Zeeman_splitting\');
addpath(genpath(zeeman_folder));

% Physical constants
const = init_constants();

% Parameters
opt.mg_range = 2;%0:2;
opt.e_J_range = [2,3];
B_range = linspace(0,20,50);
g_level = '2_3P_2';
g_state = '2_3P_2_2';

e_term = '5_3D';
e_manifold  = e_term(1:3);
% Generate the lookup table
% lines = z_level2term(B,g_level,e_term,const);
lines = z_level2manifold(B_range,g_level,e_manifold,const);

sfigure(123);
clf;
plot_level2term(lines,g_level,e_term,const,opt)
hold on
% errorbar(B_exp'.*ones(size(f_exp)),1e6*f_exp,1e6*(w_exp+f_err),'kx')
errorbar(B_exp'.*ones(size(f_single)),1e6*f_single,1e6*(w_single),'kx')
% plot(B_exp,1e6*f_add,'kx')

xlabel('Magnetic field strength (G)')
ylabel('Frequency (Hz)')

