%% Data analysis for Helium spectroscopy

clear all;
% Remove old data dirs from path


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% GETTING STARTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_dir = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
e_state = '5^3D_1';
this_folder = fileparts(which(mfilename));
core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
zeeman_folder = fullfile(fileparts(this_folder),'Zeeman_splitting\');
addpath(genpath(this_folder));
addpath(genpath(core_folder));
addpath(genpath(zeeman_folder));
addpath(data_dir)
fwtext('')
fwtext('STARTING ANALYSIS')
fwtext('')

header({0,'Setting up configs...'})
% Declare useful constants
hebec_constants
% initialize variables
% Note that master config also calls a local override if it exists
opts = master_transition_config(data_dir);
header({1,'Done.'})

% lopt_file = fullfile(data_dir,'local_opts.m');
% if exist(lopt_file,'file')==2
%     lopts = function_handle(lopt_file);
%     opts = lopts(opts,data_dir);
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Import the analog log files
data.ai = ai_log_import(opts);
%% Import LabView log
data.lv = import_lv_log(opts);
% %% Import wavemeter logs
data.wm = wm_log_import(opts);
% %% Import TDC files
data.tdc = import_mcp_tdc(opts);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Match the timestamps    
data.sync = match_timestamps(data,opts);

%% Create a calibration model
data.cal = make_calibration_model(data,opts);

%% Mask out shots which failed
data.check = check_for_errors(data,opts);

%% Break data into categories
data.cat = categorize_shots(data,opts);

%% Grouping by wavelength 
data = bin_by_wavelength(data,opts);

%% Peak detection
data = auto_peak_detect(data,opts);

%% Fit the detected peaks
data = fit_detected_peaks(data,opts);

%% Zeeman shift correction
% data = zeeman_correction(data,opts);
% zeeman_structured
% hold on
% %
% for cidx=1:2
%     freqs = data.cat{cidx}.pfits.lor_prms(:,1)*1e6;
%     Bval = opts.Bfield(cidx)*ones(size(freqs));
%     plot(Bval,freqs,'kx')
% end
    

%% Presentation plots


%% Save to output
header({0,'Saving output...'})
out_data.data = data.cat;
out_data.options = opts;
filename = fullfile(opts.out_dir,'output_and_options.m');
save(filename,'out_data','-v7.3')
fwtext('All Done!')
