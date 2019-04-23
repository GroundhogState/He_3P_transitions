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

clear all;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% GETTING STARTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add all subfolders to the path

data_dir = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190417_5^1D2_mj_1_itc_both_qwp_146_overnight\';

this_folder = fileparts(which(mfilename));
core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
zeeman_folder = fullfile(fileparts(this_folder),'Zeeman_\');
addpath(genpath(this_folder));
addpath(genpath(core_folder));
addpath(genpath(zeeman_folder));
fwtext('')
fwtext('STARTING ANALYSIS')
fwtext('')

%% Setting up
header({0,'Setting up configs...'})
% Declare useful constants
hebec_constants
% initialize variables

opts = master_transition_config();
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
tdc_time = data.tdc.time_create_write(:,2);
lv_time = data.lv.time;
wm_time = data.wm.feedback.posix_time;

all_setpts = unique(data.lv.setpoint);
lv_set_range = [min(all_setpts)/1e6,max(all_setpts)/1e6]; %MHz
mid_setpt = mean(lv_set_range); %MHz

data.start_time = min(lv_time);
num_files = length(lv_time);

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


%% Grouping by wavelength 
data = bin_by_wavelength(data,opts);
%% Peak detection
data = auto_peak_detect(data,opts);

%% Fitting goes here

%% Fit the detected peaks
data = fit_detected_peaks(data,opts);

%% Presentation plot



%% Try looking for hidden peaks?

%% Zeeman analysis

%%
% figure(1)
% clf
% zeeman_structured()
% f_offset = 3e6;
% B_vals = [17.8,10];
% B_err = [1.5,1];
% % B_vals = [17.8,11.951];
% %B_vals = [14.5,8.65];
% peak_idxs = [1,2,3,4,5;
%             1,3,4,4,5];
% % peak_idxs = [1,2,3,4,5;
% %             1,3,4,4,5];
% X = [];
% Y = [];
% Y_err = [];
% X_err = [];
% for cat_idx = 1:numel(data.cat)
% %     X = [X,B_vals(cat_idx)*ones(size(data.cat{cat_idx}.peaks.freqs))];
% %     X_err = [X_err,B_err(cat_idx)*[1,1,1,1,1]'];
% %     Y = [Y,data.cat{cat_idx}.peaks.freqs*1e6+f_offset];
% %     Y_err = [Y_err,1e6*[5,10,5,5,5]'];
%     X = B_vals(cat_idx)*ones(size(data.cat{cat_idx}.peaks.freqs));
%     X_err = [X_err,B_err(cat_idx)*[1,1,1,1,1]'];
%     Y = data.cat{cat_idx}.peaks.freqs*1e6+f_offset;
%     plot(X',Y','kx')
%     hold on
% end
% Y = [
%      154.46,176.78
%      187.50,187.27
%      187.50,201.50
%      216.99,218.98
%     ];
% Y = (Y + 744396000).*1e6+f_offset;
% X = ones(size(Y,1),2).*B_vals;
% plot(X',Y','k*')
% Y = [189.73,201.75
%      218.77,223.78
%      233.98,223.78
%     ];
% Y = (Y + 744396000).*1e6+f_offset;
% X = ones(size(Y,1),2).*B_vals;
% plot(X',Y','r+')
% title('Combining theory and data')

%% Fit the peaks
% data.fits = fit_found_peaks(data,opts);
%% Save to output
header({0,'Saving output...'})
out_data.data = data.cat;
out_data.options = opts;
filename = fullfile(opts.out_dir,'output_and_options.m');
save(filename,'out_data','-v7.3')
fwtext('All Done!')

%% Plotting

%   
% % % Compute things
% plot_pred = opts.pred_freq.*1e3-data.fits.gaus_fit(1);
% 
% plot_set_X = 2*sort(all_setpts);
% plot_set_X = plot_set_X(~isnan(plot_set_X))/1e6; %MHz
% 
% plot_fit_X = linspace(min([data.fits.freqs_centred,plot_pred-1]),max([data.fits.freqs_centred,plot_pred+1]),1e4);
% plot_fit_Y_lor = data.fits.lor_mdlfun(data.fits.lor_cen_prms,plot_fit_X);
% plot_fit_Y_gaus = data.fits.gaus_mdlfun(data.fits.gaus_cen_prms,plot_fit_X);
% 
% plot_sig_X = data.freq_stats.freq-data.fits.lorz_fit(1);
% plot_sig_Y = data.freq_stats.sig_cal;
% plot_sig_Y_err = data.freq_stats.sig_err;
% plot_sig_Y_err = plot_sig_Y_err.*~isnan(plot_sig_Y_err);
% 
% % Plotting the result
% f3=sfigure(501);
% clf;
% errorbar(plot_sig_X,plot_sig_Y,plot_sig_Y_err,'.')
% hold on
% plot(plot_fit_X,plot_fit_Y_lor,'r-')
% plot(plot_fit_X,plot_fit_Y_gaus,'k--')
% plot(plot_pred.*[1,1],[0,2],'b-','LineWidth',2.0)
% xlim(1.1*[min([plot_sig_X;plot_pred-1]),max([plot_sig_X;plot_pred+1])])
% ylim([0,1.15])
% xlabel(sprintf('f - %6f (MHz)',data.fits.lorz_fit(1)))
% ylabel('Atom Number Ratio')
% legend('Frequency-binned data','Lorentzian fit','Gaussian fit','Theory Value','Location','Best')
% 
% suptitle([opts.tr_name,sprintf(' absorption peak @ %.2f±(%.3f) MHz, FWHM~%.2fMHz',...
%     data.fits.lorz_fit(1),data.fits.lorz_fit(2),data.fits.lorz_fit(3))])
% title(sprintf('Vacuum wavelength %.10f nm ± %.6f fm, 2*%u shots',...
%     data.fits.lorz_fit(5),1e6*data.fits.lorz_fit(6),length(data.sync.msr.tdc_time)))
% 
% % Saving 
% filename3 = fullfile(opts.out_dir,sprintf('%s_spectrum',opts.tr_name));
% saveas(f3,[filename3,'.fig']);
% saveas(f3,[filename3,'.png'])
% 
% 
% fprintf('All done!\n')
% header({1,'Done.'})

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% WRITE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fwtext(' RESULTS ')
% header({1,'Peak data'})
% fprintf('Transition Frequency             %.2f±(%.3f) MHz\n',data.fits.gaus_fit(1),data.gaus_fit(2))
% fprintf('Peak width                       %.3f±(%.3f) MHz\n',data.gaus_fit(3),data.gaus_fit(4))
% fprintf('Transition Wavelength            %.6f±(%.7f) nm\n',data.gaus_fit(5),data.gaus_fit(6))
% header({1,'Lorentzian fit'})
% fprintf('Transition Frequency             %.2f±(%.3f) MHz\n',data.fits.lorz_fit(1),data.lorz_fit(2))
% fprintf('Peak width                       %.3f±(%.3f) MHz\n',data.lorz_fit(3),data.lorz_fit(4))
% fprintf('Transition Wavelength            %.6f±(%.7f) nm\n',data.lorz_fit(5),data.lorz_fit(6))
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Commented out for now; eventually the fit parameters & uncerts can be written out also
%% Settings

