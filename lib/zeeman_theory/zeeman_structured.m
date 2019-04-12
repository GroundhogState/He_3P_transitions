% Note that J ranges from |L-S| to L+S
% A combination of S and L is called a term
% S, L, J specify a level
% S, L, J, and m_J specify a state

%% OK, this is pretty functional!
% Take a break to get the analysis working. Then revisit this and consider what you want to pull from it.

%% init
disp('=== Zeeman splittings ===')
% clear all;
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
    % Updated c.f. Drake's email
const.f_table.g_2_3P_2.e_5_3D_3 = 744396208.36e6;%744.39620968e12;
const.f_table.g_2_3P_2.e_5_3D_2 = 744396227.58e6;% 744.39622889e12;
const.f_table.g_2_3P_2.e_5_3D_1 = 744396511.14e6 ;%744.39651246e12; 
% Singlet-triplet transitions
const.f_table.g_2_3P_2.e_5_1S_0 = 1e9*const.c/406.8886971706; % Can't observe from our pump state :( 
const.f_table.g_2_3P_2.e_5_1P_1 = 1e9*const.c/402.322271224483;  % Should be visible with sigma-
const.f_table.g_2_3P_2.e_5_1D_2 = 744430343.14e6;% Spotted; should be able to get two lines


%Fitted values for the 5^3D's
const.f_table.g_2_3P_2.e_5_3D_3 = 744.39620836e12;
const.f_table.g_2_3P_2.e_5_3D_2 = 744.39622758e12;
const.f_table.g_2_3P_2.e_5_3D_1 = 744.39651114e12;



opt.mg_range = 2;%0:2;

% Parameters
B = linspace(0,30,50);
g_level = '2_3P_2';
g_state = '2_3P_2_2';
e_term = '5_1D_2';
e_manifold  = e_term(1:3);
% Generate the lookup table
lines = z_level2manifold(B,g_level,e_manifold,const);

plot_level2term(lines,g_level,e_term,const,opt)
% plot_state2manifold(lines,g_state,e_manifold,const,opt)
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
        const.eps = 8.854*10^(-12);
        const.e = 1.6021e-19;
        % Notation & lookup
        const.terms = {'S','P','D','F','G'};

end








