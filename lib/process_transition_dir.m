function out_data = analyse_transition_dir(data_dir)

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

   

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% GETTING STARTED
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add all subfolders to the path
    this_folder = fileparts(which(mfilename));
    core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
    addpath(genpath(this_folder));
    addpath(genpath(core_folder));
    fwtext('')
    fwtext('STARTING ANALYSIS')
    fwtext('')

    %% Setting up
    header({0,'Setting up configs...'})
    % Declare useful constants
    hebec_constants
    % initialize variables
    if ~strcmp(data_dir(end),'\'), data_dir = [data_dir,'\']; end
    opts = master_transition_config();
    % Configuring I/O directories & setting local configs
    opts = dir_config(opts,data_dir);
    
    header({1,'Done.'})

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTING DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Import the analog log files
    data.ai = ai_log_import(opts);
    %% Import LabView log
    data.lv = import_lv_log(opts);
    %% Import wavemeter logs
    data.wm = wm_log_import(opts);
    %% Import TDC files
    data.tdc = import_mcp_tdc(opts);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Declare global variables
    fwtext('PROCESSING IMPORTS')
    data.start_time = min(data.lv.time);

    %% %%%%%%%%% COMPUTATION ON THE DATA

    %% Match the timestamps    
    data.sync = match_timestamps(data,opts);
    % 
    % %% Ignore shots with errors
    % data.check = check_error_logs(data,opts);

    %% Create a calibration model
    data.sync.msr.calib = make_calibration_model(data,opts);

    %% Break data into categories
    data.cat = categorize_shots(data,opts);

    %% Peak detection
    data = auto_peak_detect(data,opts);

    %% Fitting goes here

    %% Fit the detected peaks
    data = fit_detected_peaks(data,opts);

    %% Grouping by wavelength 
    data = bin_by_wavelength(data,opts);

    %% Try looking for hidden peaks

    %% Zeeman analysis


    %% Fit the peaks
    % data.fits = fit_found_peaks(data,opts);
    %% Save to output
    header({0,'Saving output...'})
    out_data.data = data.cat;
    out_data.options = opts;
    filename = fullfile(opts.out_dir,'output_and_options.m');
    save(filename,'out_data','-v7.3')
    header({1,'Done!'})


end


