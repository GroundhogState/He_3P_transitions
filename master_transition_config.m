function [opts,const] = master_transition_config()

%% Frequently adjusted quantities

% These variables are set here for quick access to quantities passed to
% functions in fields defined later in this code
% opts.ignorefiles = 1028:1305;
peak_cutoff_thresh = 0.1;
peak_smooth_width = 15;
peak_saturation_threshold = 0.975;


% Experimental parameters
opts.probe_set_pt=0.4;
opts.trig_ai_in=20;
opts.trig_dld=20.3;
opts.dld_aquire=4;
opts.aquire_time=4;
opts.aom_freq=189;%190*1e6;%Hz %set to zero for comparison with previous data runs
% Experimental constants


%% LabView import



opts.lv.plots = true;

%% Analog import
opts.ai.num_files = nan;

% opts.ai.force_recalc_import = true;
opts.ai.log_name='log_analog_in_';
opts.ai.verbose = 1;
opts.ai.plots=1;
opts.ai.post_fun = @transition_post_ai;
% Options for global import
opts.ai.cache_import.verbose=0;
opts.ai.cache_import.force_cache_load=false;
opts.ai.cache_import.force_recalc=false;
% Options for single imports (unnecessary unless they're cached, ill-advised usually)
opts.ai.cache_single_import = false;
opts.ai.cache_single.verbose=0;
opts.ai.cache_single.force_recalc=false;
opts.ai.cache_single.path_directions={1,'dir'};
opts.ai.args_single.cmp_multiplier_disp=50; %multiplier to display the compressed data better





%% Wavemeter log importing

opts.wm.force_reimport=false;
opts.wm.num_logs = nan;
opts.wm.plots = true;
opts.wm.plot_all=1;
opts.wm.plot_failed=false;


opts.wm.wm_log_name='log_wm_';

opts.wm.cache_import.verbose=0;
opts.wm.cache_import.force_recalc=0;

opts.wm.cache_import.save_compressed=true;%needed otherwise save takes a very long time
opts.wm.cache_import.path_directions={1,'dir'};

opts.wm.plot_all=true;
opts.wm.plot_failed=false;
opts.wm.force_reimport=false;



%% TDC import


opts.tdc.plots = true;
opts.tdc.file_name='d';
opts.tdc.force_load_save=false;   %takes precidence over force_reimport
opts.tdc.force_reimport=false;
opts.tdc.force_forc=false;
opts.tdc.dld_xy_rot=0.61;
%Should probably try optimizing these
tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
opts.tdc.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.max_runtime=inf;%inf%cut off the data run after some number of hours, should bin included as its own logic not applied to the atom number ok


%% Error handling
opts.check.wm_tolerance = 2;
opts.check.num_cal_bins = 20;

%% Peak detection
opts.peak.plot = true;
opts.peak.cutoff_thresh = .15;
opts.peak.smooth_width = 15;
opts.peak.saturation_threshold = 0.975;


%% Plotting


%% Physical constants

const.mu = 9.27e-28; %J/G
const.h = 6.63e-34;
const.hbar = const.h/(2*pi);
const.f_mu = const.mu/const.h;
const.w_mu = const.mu/const.hbar;
const.c = 299792458;
% Notation & lookup
const.terms = {'S','P','D','F','G'};
%% Reference values
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

%Fitted valuse for the 5^3D's
const.f_table.g_2_3P_2.e_5_3D_3 = 744.39620836e12;
const.f_table.g_2_3P_2.e_5_3D_2 = 744.39622758e12;
const.f_table.g_2_3P_2.e_5_3D_1 = 744.39651114e12;


%% Experimental constants
const.fall_time=0.417;
const.global.qe=0.09;

%% Append to options
opts.const = const;


end



