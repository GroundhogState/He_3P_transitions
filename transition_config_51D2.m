function opts = transition_config()

opts.dir = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190401_singlet_5^1D2_fine_scan\';
opts.dir = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190405_5^1D_2_cal_qwp_236_stage_1_ITC\';
% opts.dir = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190405_5^1D_2_cal_qwp_146_stage_1_ITC\'
opts.transition = '5^1D_2';
opts.tr.tr_name = opts.transition;

opts.probe_set_pt=0.4;
opts.ritz_wl = 402.7151538; % NIST Ritz wavelength, nm
opts.ritz_freq = 299792458/opts.ritz_wl; %GHz
opts.obs_wl = 412.1978; %NIST observed vacuum wl, nm
opts.obs_freq = 299792458/opts.obs_wl; %GHz

opts.pred_wl = opts.ritz_wl;
opts.pred_freq = opts.ritz_freq+2.29;
%opts.pred_freq = 744.396188e3;

opts.global.fall_time=0.417;
opts.global.qe=0.09;

opts.wm.plot_all=1;

opts.wm.plot_failed=false;
opts.trig_ai_in=20;
opts.trig_dld=20.3;
opts.dld_aquire=4;
opts.aquire_time=4;

%% Generate out directory
%set up an output dir %https://gist.github.com/ferryzhou/2269380
if (exist([opts.dir,'out'], 'dir') == 0), mkdir([opts.dir,'out']); end
%make a subfolder with the ISO timestamp for that date
opts.out_dir=sprintf('%sout\\%s\\',...
    opts.dir,datestr(datetime('now'),'yyyymmddTHHMMSS'));
if (exist(opts.out_dir, 'dir') == 0), mkdir(opts.out_dir); end


%% Analog import
opts.ai.num_files = nan;
opts.ai.dir=opts.dir;
opts.ai.out_dir = opts.out_dir;
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
opts.ai.cache_single.mock_working_dir=opts.dir;
opts.ai.cache_single.path_directions={1,'dir'};
opts.ai.args_single.cmp_multiplier_disp=50; %multiplier to display the compressed data better


opts.ai.pd.diff_thresh=0.1;
opts.ai.pd.std_thresh=0.1;
opts.ai.pd.time_start=0.2;
opts.ai.pd.time_stop=2;
opts.ai.sfp.num_checks=20; %how many places to check that the laser is single mode
opts.ai.sfp.thresh_cmp_peak=20e-3; %theshold on the compressed signal to be considered a peak
opts.ai.sfp.peak_dist_min_pass=4.5;%minimum (min difference)between peaks for the laser to be considered single mode
opts.ai.plot.all=false;
opts.ai.plot.failed=false;
opts.ai.time_match_valid=8; %how close the predicted start of the shot is to the actual
opts.ai.scan_time=1/20; %fast setting 1/100hz %estimate of the sfp scan time,used to set the window and the smoothing
%because im only passing the ai feild to aviod conflicts forcing a reimport i need to coppy these feilds
opts.ai.trig_ai_in=opts.trig_ai_in;
opts.ai.trig_dld=opts.trig_dld;
opts.ai.dld_aquire=opts.dld_aquire;
opts.ai.aquire_time=opts.dld_aquire;
opts.ai.tdc_override=1;


%% Wavemeter log
opts.wm.dir=opts.dir;
opts.wm.out_dir = opts.out_dir;
opts.wm.force_reimport=false;
opts.wm.num_logs = nan;
opts.wm.plots = true;

wm_log_name='log_wm_';
wm_logs=dir([opts.wm.dir,wm_log_name,'*.txt']);
opts.wm.names={wm_logs.name};

opts.wm.cache_import.verbose=0;
opts.wm.cache_import.force_recalc=0;
opts.wm.cache_import.mock_working_dir=opts.dir;
opts.wm.cache_import.save_compressed=true;%needed otherwise save takes a very long time
opts.wm.cache_import.path_directions={1,'dir'};

opts.wm.plot_all=true;
opts.wm.plot_failed=false;
opts.wm.force_reimport=false;

opts.wm.time_pd_padding=4; %check this many s each side of probe
opts.wm.time_blue_padding=1; %check this many seconde each side of probe
opts.wm.time_probe=3;
opts.wm.ecd_volt_thresh=0.5;

opts.wm.red_sd_thresh=50; %allowable standard deviation in MHz
opts.wm.red_range_thresh=50; %allowable range deviation in MHz
opts.wm.rvb_thresh=20; %allowable value of abs(2*red-blue)

%% LV import

opts.lv.dir = opts.dir;
opts.lv.out_dir = opts.out_dir;
opts.lv.plots = true;


%% TDC import

opts.tdc.dir = opts.dir;
opts.tdc.out_dir = opts.out_dir;
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

opts.aom_freq=0;%190*1e6;%Hz %set to zero for comparison with previous data runs


%% Plotting
opts.tr.plot = 1; 
opts.tr.out_dir = opts.out_dir;
opts.tr.pred_freq = opts.pred_freq;
opts.tr.pred_wl = opts.pred_wl;
opts.tr.num_cal_bins = 50;
opts.tr.wm_tolerance = 10; %MHz
opts.tr.num_freq_bins = 30; %MHz
opts.tr.freq_bin_size = 0.25;
opts.tr.aom_freq = 189; %MHz, blue


end