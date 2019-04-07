% Note that J ranges from |L-S| to L+S
% A combination of S and L is called a term
% S, L, J specify a level
% S, L, J, and m_J specify a state

%% OK, this is pretty functional!
% Take a break to get the analysis working. Then revisit this and consider what you want to pull from it.

%% init
disp('=== Zeeman splittings ===')
clear all;
this_folder = fileparts(which(mfilename));
core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
addpath(genpath(this_folder));
addpath(genpath(core_folder));


% Physical constants
const = init_constants();
% Table of field-free freqs defined here for readability
% Read N_(2S+1)L_J 
% The 412nm benchmark
const.f_table.g_2_3P_2.e_5_3S_1 = 727.3032446e12;
% Misc transitions - what do the stars mean?
const.f_table.g_2_3P_2.e_5_3P_0 = 1e9*const.c/404.628937550957;
const.f_table.g_2_3P_2.e_5_3P_1 = 1e9*const.c/404.629844755577;
const.f_table.g_2_3P_2.e_5_3P_2 = 1e9*const.c/404.629918705477; 
% Historically controversial transitions
const.f_table.g_2_3P_2.e_5_3D_3 = 744.39620968e12;
const.f_table.g_2_3P_2.e_5_3D_2 = 744.39622889e12;
const.f_table.g_2_3P_2.e_5_3D_1 = 744.39651246e12; 
% Singlet-triplet transitions
const.f_table.g_2_3P_2.e_5_1S_0 = 1e9*const.c/406.8886971706;
const.f_table.g_2_3P_2.e_5_1P_1 = 1e9*const.c/402.322271224483;
const.f_table.g_2_3P_2.e_5_1D_2 = 744.43034335e12; % 402.7nm

opt.prm = [];

% Parameters
B = linspace(0,25,100);
g_level = '2_3P_2';
g_state = '2_3P_2_2';
e_term = '5_3D';
% Generate the lookup table
% lines = z_level2term(B,g_level,e_term,const);

% Plot things, looping over upper levels and then types of transition
% function plot_state2manifold


figure(1)
clf
    % kind of built exclusively for each other atm
    e_manifold = '5_3'; % singlet manifold
    lines = z_level2manifold(B,g_level,e_manifold,const);
    plot_state2manifold(lines,g_state,e_manifold,const,opt)
% plot_level2term(lines,g_level,e_term,const,opt)
% plot_state2term(lines,g_state,e_term,const,opt)
xlabel('Magnetic field strength (G)')
ylabel('Frequency (Hz)')
% g_label = [strrep(g_state(1:4),'_','^'),g_state(5:6),strrep(g_state(7:end),'_',',m_j=')];
% title(sprintf('Zeeman lines of the %s -> %s transition',g_label,strrep(e_term,'_','^')))

disp('  = Done =  ')


%% functions
% Could add a term2term function as well, but we're pretty sure we just
% live in the 2^3P_2_mj_2...

function const = init_constants()

        const.mu = 9.27e-28; %J/G
        const.h = 6.63e-34;
        const.hbar = const.h/(2*pi);
        const.f_mu = const.mu/const.h;
        const.w_mu = const.mu/const.hbar;
        const.c = 299792458;
        % Notation & lookup
        const.terms = {'S','P','D','F','G'};

end








