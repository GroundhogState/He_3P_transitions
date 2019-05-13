%% Making a level diagram
clear all
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
d_opts = master_transition_config('null');
const = opts.const;
cooling_wl = 1083.331e-9; %nm
cool_gap = const.c/cooling_wl; %Hz
ms_gap = 19.8*const.q/const.h; %Hz

e_states = {'5_3S_1','5_3D_1','5_3D_2','5_3D_3','5_1D_2'};
lvls.states = {'1_1S_0','2_3S_1','2_3P_2',e_states{:}}
nstates = numel(e_states);
lvls.gaps = zeros(nstates,1);
lvls.Ls = zeros(nstates,1);
for nstate = 1:nstates
    s_upper = e_states{nstate};
    Ls(nstate) = find(strcmp(s_upper(4),const.terms))-1;
    gaps(nstate) = const.f_table.('g_2_3P_2').(['e_',s_upper])+cool_gap+ms_gap; %Hz
end

lvls.levels = [0,ms_gap,cool_gap+ms_gap,gaps']*const.h/const.q; %eV
lvls.Ls = [0,0,1,Ls']
