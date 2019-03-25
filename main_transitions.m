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
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
fwtext('')
fwtext('STARTING ANALYSIS')
fwtext('')
data = [];
%% Setting up
header({0,'Setting up configs...'})
% Declare useful constants
hebec_constants
% initialize variables
opts = transition_config();
header({1,'Done.'})

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwtext('IMPORTING DATA')
% % IMPORT THE ANALOG INPUT
data.ai=ai_log_import(opts.ai,data);


%% PD data is in channel 3 
% Code below will develop into a post-process fn to be called on data.ai_log

%% Import LabView log
data.lv = import_lv_log(opts);
%% IMPORT WM LOG FILES
data.wm=wm_log_import(opts.wm);
% data.wm
% %% CHECK THE WM INPUTS ?
% data.wm.proc=wm_log_process(opts,data);
% clear('sub_data')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwtext('PROCESSING IMPORTS')
header({1,'Processing...'})
data.tr = transition_process(data,opts.tr);
% transition_plots(data,opts.plt)

fprintf('All done!\n')

header({1,'Done.'})
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PRESENTING RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fwtext('ALL DONE!!!')
%% Settings
