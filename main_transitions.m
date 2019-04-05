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
opts = transition_config_51D2();
header({1,'Done.'})

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import the analog log files
%data.ai = ai_log_import(opts.ai);
%% Import LabView log
data.lv = import_lv_log(opts.lv);
%% Import wavemeter logs
data.wm = wm_log_import(opts.wm);
%% Import TDC files
data.tdc = import_mcp_tdc(opts.tdc);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwtext('PROCESSING IMPORTS')
header({1,'Processing...'})
% data.tr = transition_process(data,opts.tr);
% data = bec_transition_process(data,opts.tr);


%% %%%%%%%%%%%% PRECONDITIONING THE DATA 
%     num_files = numel(data.tdc.shot_num);
tdc_time = data.tdc.time_create_write(:,2);
lv_time = data.lv.time;
wm_time = data.wm.feedback.posix_time;

all_setpts = unique(data.lv.setpoint);
lv_set_range = [min(all_setpts)/1e6,max(all_setpts)/1e6]; %MHz
mid_setpt = mean(lv_set_range); %MHz

data.tr.start_time = min(lv_time);
num_files = length(lv_time);

%% %%%%%%%%% COMPUTATION ON THE DATA

%% Match the timestamps    
data.tr.sync = match_timestamps(data,opts.tr);

%% Create a simple calibration model
data.tr.sync.msr.calib = make_calibration_model(data,opts.tr);

%% Grouping by wavelength (Independent variable)
data.tr.freq_stats = bin_by_wavelength(data,opts.tr);


%% Fitting the data

% Finding the peaks
opts_tr.cutoff_thresh = 0.1;
opts_tr.smooth_width = 20;
opts_tr.saturation_threshold = 0.975;
data.tr.peaks = find_spectral_peaks(data.tr,opts_tr);

% Fit the peaks
data.tr.fits = fit_peaks(data,opts_tr);

%% Plotting
if opts_tr.plot
  
    % % Compute things
    plot_pred = opts_tr.pred_freq.*1e3-data.tr.fits.gaus_fit(1);
    
    plot_set_X = 2*sort(all_setpts);
    plot_set_X = plot_set_X(~isnan(plot_set_X))/1e6; %MHz

    plot_fit_X = linspace(min([data.tr.fits.freqs_centred,plot_pred-1]),max([data.tr.fits.freqs_centred,plot_pred+1]),1e4);
    plot_fit_Y_lor = data.tr.fits.lor_mdlfun(data.tr.fits.lor_cen_prms,plot_fit_X);
    plot_fit_Y_gaus = data.tr.fits.gaus_mdlfun(data.tr.fits.gaus_cen_prms,plot_fit_X);
        
    plot_sig_X = data.tr.freq_stats.freq-data.tr.fits.lorz_fit(1);
    plot_sig_Y = data.tr.freq_stats.sig_cal;
    plot_sig_Y_err = data.tr.freq_stats.sig_err;
    plot_sig_Y_err = plot_sig_Y_err.*~isnan(plot_sig_Y_err);
    
    % Plotting the result
    f3=sfigure(501);
    clf;
    errorbar(plot_sig_X,plot_sig_Y,plot_sig_Y_err,'.')
    hold on
    plot(plot_fit_X,plot_fit_Y_lor,'r-')
    plot(plot_fit_X,plot_fit_Y_gaus,'k--')
    plot(plot_pred.*[1,1],[0,2],'b-','LineWidth',2.0)
    xlim(1.1*[min([plot_sig_X;plot_pred-1]),max([plot_sig_X;plot_pred+1])])
    ylim([0,1.15])
    xlabel(sprintf('f - %6f (MHz)',data.tr.fits.lorz_fit(1)))
    ylabel('Atom Number Ratio')
    legend('Frequency-binned data','Lorentzian fit','Gaussian fit','Theory Value','Location','Best')

    suptitle([opts_tr.tr_name,sprintf(' absorption peak @ %.2f±(%.3f) MHz, FWHM~%.2fMHz',...
        data.tr.fits.lorz_fit(1),data.tr.fits.lorz_fit(2),data.tr.fits.lorz_fit(3))])
    title(sprintf('Vacuum wavelength %.10f nm ± %.6f fm, 2*%u shots',...
        data.tr.fits.lorz_fit(5),1e6*data.tr.fits.lorz_fit(6),length(data.tr.sync.msr.tdc_time)))
    
    % Saving 
    filename3 = fullfile(opts_tr.out_dir,sprintf('%s_spectrum',opts_tr.tr_name));
    saveas(f3,[filename3,'.fig']);
    saveas(f3,[filename3,'.png'])
end

fprintf('All done!\n')
header({1,'Done.'})

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% WRITE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwtext(' RESULTS ')
header({1,'Gaussian fit'})
fprintf('Transition Frequency             %.2f±(%.3f) MHz\n',data.tr.fits.gaus_fit(1),data.tr.gaus_fit(2))
fprintf('Peak width                       %.3f±(%.3f) MHz\n',data.tr.gaus_fit(3),data.tr.gaus_fit(4))
fprintf('Transition Wavelength            %.6f±(%.7f) nm\n',data.tr.gaus_fit(5),data.tr.gaus_fit(6))
header({1,'Lorentzian fit'})
fprintf('Transition Frequency             %.2f±(%.3f) MHz\n',data.tr.fits.lorz_fit(1),data.tr.lorz_fit(2))
fprintf('Peak width                       %.3f±(%.3f) MHz\n',data.tr.lorz_fit(3),data.tr.lorz_fit(4))
fprintf('Transition Wavelength            %.6f±(%.7f) nm\n',data.tr.lorz_fit(5),data.tr.lorz_fit(6))
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Commented out for now; eventually the fit parameters & uncerts can be written out also
% header({0,'Saving output'})
% out_struct.data = data;
% out_struct.opts = opts;
% 
% % Trim the raw data
% out_struct.data.tdc = rmfield(out_struct.data.tdc,'counts_txy');
% % WM and AI logs would be huge
% out_struct.data = rmfield(out_struct.data,'wm');
% out_struct.data = rmfield(out_struct.data,'ai');
% outfilename = fullfile(opts.out_dir,['results_',datestr(datetime('now'),'yyyymmddTHHMMSS'),'.mat']);
% save(outfilename,'out_struct','-v7.3')
% header({0,'Saved'})
% fwtext('ALL DONE!!!')
%% Settings

