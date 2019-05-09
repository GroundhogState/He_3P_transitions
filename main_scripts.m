%% Data analysis for Helium spectroscopy

% clear all;
% Remove old data dirs from path

data = [];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% GETTING STARTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% add all subfolders to the path
% % Setting up
data_dir = '\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\2019_He_Spectroscopy\spectroscopy_data\final_scans\20190429_5^3D_2_3_qwp_142\';
% data_dir = 'Z:\EXPERIMENT-DATA\2019_He_transitions\20190409_5^3D2_3D3_qwp_146_two_stage\';


this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
addpath(genpath(core_folder));
zeeman_folder = fullfile(fileparts(this_folder),'Zeeman_splitting\');
addpath(genpath(zeeman_folder));

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

% lopt_file = fullfile(data_dir,'local_opts.m');
% if exist(lopt_file,'file')==2
%     opts = local_opts(opts);
% end

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
%% Match the timestamps    
data.sync = match_timestamps(data,opts);

% %% Create a calibration model
data.cal = make_calibration_model(data,opts);

% Mask out shots which failed
data.check = check_for_errors(data,opts);

%% Break data into categories
data.cat = categorize_shots(data,opts);
 
% %% Peak detection
data = auto_peak_detect(data,opts);

%% Fit the detected peaks
data = fit_detected_peaks(data,opts);
% data = fancy_fits(data,opts);
% 
%% Zeeman shift correction
data = zeeman_correction(data,opts);

% %% Presentation plots

lfun = @(p,x) p(1)*p(3)./((x-p(2)).^2 + (p(3)).^2);
cli_header({0,'Forming final plots...'});
ncat = numel(data.cat);


e_level = opts.e_state(1:6);
fmt_name = strrep(e_level,'^','_');
f_pred = opts.const.f_table.g_2_3P_2.(sprintf('e_%s',fmt_name))/1e6; %MHz

fprintf('Predicted vacuum freq %.3f MHz\n',f_pred)
if ~iscell(opts.e_state)       
    num_pks = 1;
else
    num_pks = numel(opts.e_state);
end

sfigure(7777);
clf;
for cidx =1:ncat
   for pidx=1:num_pks
       cdata = data.cat{cidx};
       B = opts.Bfield(cidx); 
       f_cen = cdata.zeeman.corrected(pidx);

       fnum = 7777+cidx;
       l_mdl = cdata.pfits.lorz{1};
       X=cdata.spec.freq-cdata.pfits.lor_prms(1);
       Xfit = linspace(min(X),max(X),100);
       Yfit = lfun(l_mdl.Coefficients.Estimate,Xfit);

       [ysamp,yci]=predict(l_mdl,Xfit','Prediction','observation'); 

       freq_stat_err = cdata.zeeman.shift_unc(pidx) + cdata.pfits.lor_err(pidx);
                    % Z shift unc + fit err + power dep? +  

       % to do: Include error from zeeman shift
       fprintf('Peak %u fitted value:   %.2f(%.2f) MHz\n',cidx,cdata.zeeman.corrected(pidx),freq_stat_err)
       fprintf('         theory diff:   %.2f MHz\n',cdata.zeeman.corrected(pidx)-f_pred)

       f1=sfigure(fnum);
       clf;
       plot([f_pred-f_cen,f_pred-f_cen],[-1e4,0],'Color',[0.8500 0.3250 0.0980])
       hold on
       plot(X,cdata.spec.signal,'.','MarkerSize',15,'Color',[0.8350 0.1780 0.1840])
       plot(Xfit,yci,'k:','LineWidth',1.)
       plot(Xfit,Yfit,'b')

       ylabel('N loss')
       xlabel(sprintf('Frequency-%.2f (MHz)',cdata.pfits.lor_prms(pidx,1)))
       title(sprintf('%s transition, |B|=%.2f, \\Delta=%.2f MHz',opts.e_state,B,cdata.zeeman.shift(pidx)))
       legend({'Theoretical value'},'location','NorthWest')
       set(gca,'FontSize',10,'FontName','cmr10')
       set(gcf,'color','w');
       imname = sprintf('plot_pretty_%u',cidx);
       filename1 = fullfile(opts.out_dir,imname);
       saveas(f1,[filename1,'.fig']);
       saveas(f1,[filename1,'.png'])

       sfigure(7777);
       plot(cdata.spec.freq-cdata.zeeman.shift(pidx),cdata.spec.signal,'k*')
       hold on
       plot([f_pred,f_pred],[-1e4,0],'r')
   end
end
% % % RESULTS
% Transition name
% Directory
% Various stat errors
fitted_freqs = cellfun(@(x) x.zeeman.corrected(1), data.cat);
stat_err = cellfun(@(x) x.zeeman.stat_unc(1), data.cat);
fprintf('COMBINED MEASUREMENT:  %.3f(%.3f) MHz\n',mean(fitted_freqs),sum(stat_err))
fprintf('   Theory difference:  %.3f MHz\n',mean(fitted_freqs)-f_pred)

cli_header({1,'Done.'})



% %% Save to output
cli_header({0,'Saving output...'})
out_data.data = data.cat;
out_data.options = opts;
filename = fullfile(opts.out_dir,'output_and_options.m');
save(filename,'out_data','-v7.3')
fwtext('All Done!')
