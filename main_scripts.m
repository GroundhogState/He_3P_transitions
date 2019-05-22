%% Data analysis for Helium spectroscopy


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% GETTING STARTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add all subfolders to the path
% % Setting up
data_dir = '\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\2019_He_Spectroscopy\spectroscopy_data\final_scans\20190417_5^1D2_mj_1_itc_both_qwp_146_overnight\';
% data_dir = '\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\2019_He_Spectroscopy\spectroscopy_data\final_scans\20190425_5^3D1_both_ITC_1MHz_step\';
% data_dir = '\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\2019_He_Spectroscopy\spectroscopy_data\final_scans\20190429_5^3D_2_3_qwp_142\';
% data_dir = '\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\2019_He_Spectroscopy\spectroscopy_data\final_scans\20190429_5^3D_2_3_qwp_142_tight_scan_centre_peak\';
% data_dir = '\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\2019_He_Spectroscopy\spectroscopy_data\final_scans\20190429_5^3S_1_qwp_146_both_stage\';

this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
addpath(genpath(core_folder));

addpath(data_dir)
fwtext('')
fwtext('STARTING ANALYSIS')
fwtext('')

cli_header({0,'Setting up configs...'})
% Declare useful constants
hebec_constants
% initialize variables
% Note that master config also calls a local override if it exists
opts = master_transition_config(data_dir);
cli_header({1,'Done.'})

opts.ai.cache_import.force_cache_load = true;
opts.wm.cache_import.force_cache_load = true;
opts.tdc.cache_import.force_cache_load = true;

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import the analog log files
data.ai = ai_log_import(opts);
% %% Import LabView log
data.lv = import_lv_log(opts);
% %% Import wavemeter logs
data.wm = wm_log_import(opts);
%% Import TDC files
data.tdc = import_mcp_tdc(opts);
cli_header({0,'Data import complete!'})
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lopt_file = fullfile(data_dir,'local_opts.m');
if exist(lopt_file,'file')==2
    addpath(data_dir)
    fprintf('Applying local overrides...\n')
    opts = local_opts(opts,data_dir);
    rmpath(data_dir)
end
% %% Match the timestamps    
data.sync = match_timestamps(data,opts);

% %% Create a calibration model
data.cal = make_calibration_model(data,opts);

% Mask out shots which failed
data.check = check_for_errors(data,opts);

% %% Break data into categories
data.cat = categorize_shots(data,opts);
 
% %% Peak detection
data = auto_peak_detect(data,opts);

% %% Fit the detected peaks
data = fit_detected_peaks(data,opts);

% %% Zeeman shift correction
data = zeeman_correction(data,opts);

%% Presentation plots
data = present_plots(data,opts);


%% SPIT OUT RESULTS
if iscell(opts.e_state)
    num_pks=numel(opts.e_state);
else
    num_pks = 1;
end

for pidx=1:num_pks
        if ~iscell(opts.e_state)
            e_state = opts.e_state;
        else
            e_state = opts.e_state{pidx};
        end
        e_level = e_state(1:6);
        fmt_name = strrep(e_level,'^','_');
        f_pred = opts.const.f_table.g_2_3P_2.(sprintf('e_%s',fmt_name))/1e6; %MHz
    fitted_freqs = cellfun(@(x) x.zeeman.corrected(pidx), data.cat);
    stat_err = cellfun(@(x) x.zeeman.stat_unc(pidx), data.cat);
    fprintf('===COMBINED MEASUREMENT===\n')
    fprintf('%s measured value:  %.3f(%.3f) MHz\n',e_level,mean(fitted_freqs),sum(stat_err))
    fprintf('Theory difference:  %.3f MHz\n',mean(fitted_freqs)-f_pred)
end

cli_header({1,'Done.'})



% %% Save to output
cli_header({0,'Saving output...'})
out_data.data = data.cat;
out_data.options = opts;
filename = fullfile(opts.out_dir,'output_and_options.mat');
save(filename,'out_data','-v7.3')
fwtext('All Done!')
