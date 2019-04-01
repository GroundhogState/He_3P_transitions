%% Prototype data analysis for Helium spectroscopy

% Needs to:
%     * Import analog data
%     Loop over files in some dir, return struct for each of them
%         * Photodiode trace
%         * Laser set point
%         * WM log (i.e. set point)
%         * SFP trace

% TO DO
% Match calibration shots to measurement shots by timestamp (within chosen tolerance) 
% Add some fits if you like and then
% Show calibrated signal for various methods
% Implement better checks: PD set point, SFP, for instance


%   Add a field somewhere in interface log with global params
%          eg, list of setpoints, other internal configs
%             Say, on file creation, write preformatted first line
        
% Desirable: Retrieve scan parameters from interface 



%%
%     FORGOT TO CLEAR WM LOG
%     Means we'll have to sort by timestamps; would have needed to do this anyway, so good prompt
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% GETTING STARTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add all subfolders to the path

clear all
this_folder = fileparts(which(mfilename));
core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
addpath(genpath(this_folder));
addpath(genpath(core_folder));
fwtext('')
fwtext('STARTING ANALYSIS')
fwtext('')

% % Setting up
header({0,'Setting up configs...'})
% Declare useful constants
hebec_constants
% initialize variables
opts = transition_config();
header({1,'Done.'})

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import the analog log files
data.ai = ai_log_import(opts.ai);
%% Import LabView log
data.lv = import_lv_log(opts.lv);
%% Import wavemeter logs
data.wm = wm_log_import(opts.wm);
%% Import TDC files
data.tdc=import_mcp_tdc(opts.tdc);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwtext('PROCESSING IMPORTS')
header({1,'Processing...'})
% data.tr = transition_process(data,opts.tr);
data.tr = bec_transition_process(data,opts.tr);
fprintf('All done!\n')
header({1,'Done.'})
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PRESENTING RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fwtext('ALL DONE!!!')
%% Settings
