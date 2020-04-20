%% Data analysis for Helium spectroscopy

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% GETTING STARTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add all subfolders to the path
% % Setting up
% close all
clear all
hdd_dir = 'E:\Data\Spectroscopy\spectroscopy_data\final_scans';
test_dir = 'E:\Data\Spectroscopy\spectroscopy_tests\';
power_dep_dir = 'E:\Data\Spectroscopy\spectroscopy_data\probe_power_dep';

% data_dir = fullfile(hdd_dir,'\20190429_5^3S_1_qwp_146_both_stage\'); % TDC time offset -7174
% data_dir = fullfile(hdd_dir,'20190425_5^3D1_both_ITC_1MHz_step'); % TDC time offset -7176
% data_dir = fullfile(hdd_dir,'20190429_5^3D_2_3_qwp_142_tight_scan_centre_peak'); %offset -7177
% data_dir = fullfile(hdd_dir,'20190429_5^3D_2_3_qwp_142'); %offset -7177
data_dir = fullfile(hdd_dir,'\20190417_5^1D2_mj_1_itc_both_qwp_146_overnight\');% no AI !!! offset -7176
tdc_offset = 7200;

% % Power dependency checks
% data_dir = fullfile(power_dep_dir,'20190426_100mV_5^3D_1');
% data_dir = fullfile(power_dep_dir,'20190425_150mV_5^3D_1');
% data_dir = fullfile(power_dep_dir,'20190426_200mV_5^3D_1');
% data_dir = fullfile(power_dep_dir,'20190426_250mV_5^3D_1');
% tdc_offset = 0;

this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
addpath(genpath(core_folder));
addpath(genpath(fullfile(fileparts(this_folder),'Zeeman_splitting\')));

addpath(data_dir)
fwtext('')
fwtext('STARTING ANALYSIS')
fwtext('')

cli_header(0,'Setting up configs...');
% Declare useful constants
hebec_constants
% opts.const = const;
% initialize variables
% Note that master config also calls a local override if it exists
opts = master_transition_config(data_dir);
% opts.tdc_offset = tdc_offset;
cli_header(1,'Done.');

opts.ai.cache_import.force_cache_load = true;
opts.wm.cache_import.force_cache_load = true;
opts.tdc.cache_import.force_cache_load = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import the analog log files
data.ai = ai_log_import(opts);
% fprintf('===\n')
% fprintf('Avg power: %.4f +- %.4f\n',nanmean(data.ai.high_mean(data.ai.high_mean>0)),nanstd(data.ai.high_mean(data.ai.high_mean>0)))
% fprintf('peak power: %.2f +- %.2f\n',1166.7*nanmean(data.ai.high_mean(data.ai.high_mean>0)),1166.7*nanstd(data.ai.high_mean(data.ai.high_mean>0)))
% fprintf('exposure time %.3f +- %.3f\n',mean(data.ai.high_time(data.ai.high_time>0.05)),std(data.ai.high_time(data.ai.high_time>0.05))) 
%% Import LabView log
data.lv = import_lv_log(opts);
%% Import wavemeter logs
data.wm = wm_log_import(opts);
%% Import TDC files 
data.tdc = import_mcp_tdc(opts);
cli_header(0,'Data import complete!');
%% %%
% t0 = data.lv.time(1);
% n_pick = min(size(data.tdc.time_create_write,1),length(data.lv.time));
% stfig('Time diagnostics');
% clf;
% subplot(2,1,1)
% hold on
% % plot(data.ai.timestamp-t0,'b.')
% plot(data.lv.time-t0,'r.')
% plot(data.tdc.time_create_write(:,1)-t0,'k.')
% plot(data.tdc.time_create_write(:,2)-t0,'ko')
% legend('LV','Create','Write')
% subplot(2,1,2)
% plot(data.tdc.time_create_write(1:n_pick,2)-data.lv.time(1:n_pick)','k.')
% plot(data.tdc.time_create_write(1:n_pick,2)-data.lv.time(1:n_pick)','k.')
% legend('Write - LV')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lopt_file = fullfile(data_dir,'local_opts.m');
if exist(lopt_file,'file')==2
    addpath(data_dir)
    fprintf('Applying local overrides...\n');
    opts = local_opts(opts,data_dir);   
    rmpath(data_dir)
end


%% Fit temperatures
% opts.shots_chosen = 430:450;
% opts.shots_chosen = 1;
% opts.shots_chosen = 448;
% psd_opts.single_plot = false;
% psd_opts.shots_chosen = nan;
% data.psd = measure_psd(data.tdc,psd_opts);
% profile on
% data.psd = thermometry_fits(data,opts);
% data.temps = thermometry_fits(data,opts);
opts.temp.cache_opts.force_forc = false;
opts.temp.cache_opts.force_recalc = false;
opts.temp.cache_opts.force_cache_load = false;
opts.temp.shot_nums = nan;
opts.single_shot_plot = false;
data.temp = thermometry_fits(data,opts);

%%
opts.temp.plot.plt_label = '5^1D_2';
plot_thermo_fits(data,opts.temp.plot)
% profile off
% profile viewer
%% Match the timestamps
% data.sync = match_timestamps(data,opts);

% % % %% Temp analysis
% % % T_all = (abs(data.sync.shots.temp/2)*const.g0).^2 *const.mhe/const.kb;
% % % T_cal = (abs(data.sync.cal.temp/2)*const.g0).^2 *const.mhe/const.kb;
% % % lambda_db = const.h./sqrt(2*pi*const.mhe*const.kb*T_all);
% % % n0 = data.sync.shots.N_atoms.^(2/5); %at least, proportional to this
% % % 
% % % T_kelvin = (abs(data.sync.msr.temp/2)*const.g0).^2 *const.mhe/const.kb;
% % % lambda_msr = const.h./sqrt(2*pi*const.mhe*const.kb*T_kelvin);
% % % n0_msr = data.sync.msr.N_atoms; %at least, proportional to this
% % % % ah - but the number isn't a good measure any more. The condensed fraction
% % % % dropped a lot. But we have *total number* and temperature
% % % % so -> condensed fraction -> condensed number -> peak density
% % % 
% % % % PSD = peak_density * lambda^3
% % % 
% % % stfig('temp spectra');
% % % clf;
% % % subplot(2,2,1)
% % % hold on
% % % plot(T_cal(:,1),'k.')
% % % plot(T_cal(:,2),'r.')
% % % title('Calibration shots')
% % % xlabel('Shot number')
% % % ylabel('Fit temperature (K)')
% % % 
% % % subplot(2,2,2)
% % % hold on
% % % plot(data.sync.msr.probe_set,T_kelvin(:,1),'k.')
% % % plot(data.sync.msr.probe_set,T_kelvin(:,2),'r.')
% % % ylabel('Fit temperature (K)')
% % % xlabel('Frequency')
% % % title('Measurement shots')
% % % 
% % % subplot(2,2,3)
% % % hold on
% % % plot(data.sync.shots.N_atoms,T_all(:,1),'k.')
% % % plot(data.sync.shots.N_atoms,T_all(:,2),'r.')
% % % xlabel('Number of counts')
% % % ylabel('Temperature (K)')
% % % title('Number vs temperature (all shots)')
% % % 
% % % ignoremask = data.sync.msr.N_atoms < 1e4;
% % % subplot(2,2,4)
% % % hold on
% % % plot(data.sync.msr.probe_set,data.sync.msr.N_atoms,'k.')
% % % plot(data.sync.msr.probe_set(ignoremask),data.sync.msr.N_atoms(ignoremask),'rx')
% % % xlabel('Frequency')
% % % ylabel('Number')
% % % 
% % % suptitle('Full TXY fit thermometry')

% %% Create a calibration model
% data.cal = make_calibration_model(data,opts);
% 
% %% Mask out shots which failed
% data.check = check_for_errors(data,opts);
% 
% %% Break data into categories
% data.cat = categorize_shots(data,opts);
%  
% %% Peak detection
% data = auto_peak_detect(data,opts);
% 
% %% Fit the detected peaks
% data = fit_detected_peaks(data,opts);
% 
% %% Zeeman shift correction
% data = zeeman_correction(data,opts);
% 
% %% Presentation plots
% % data = present_plots(data,opts);
% 
% 
% %% SPIT OUT RESULTS
% if iscell(opts.e_state)
%     num_pks=numel(opts.e_state);
% else
%     num_pks = 1;
% end
% opts.cell_shift = -1.9;
% for pidx=1:num_pks
%         if ~iscell(opts.e_state)
%             e_state = opts.e_state;
%         else
%             e_state = opts.e_state{pidx};
%         end
%         e_level = e_state(1:6);
%         fmt_name = strrep(e_level,'^','_');
%         f_pred = opts.const.f_table.g_2_3P_2.(sprintf('e_%s',fmt_name))/1e6; %MHz
%     fitted_freqs = cellfun(@(x) x.zeeman.corrected(pidx), data.cat);
%     stat_err = cellfun(@(x) x.zeeman.stat_unc(pidx), data.cat);
%     fprintf('===COMBINED MEASUREMENT===\n')
%     fprintf('%s measured value:  %.3f(%.3f) MHz\n',e_level,mean(fitted_freqs),sum(stat_err))
%     fprintf('Theory difference:  %.3f MHz\n',mean(fitted_freqs)-f_pred)
% end
% 
% %% FINAL PLOTS
% lfun = @(p,x) p(1)*p(3)./((x-p(2)).^2 + (p(3)).^2);
% ncat = numel(data.cat);
% savedir = 'C:\Users\jaker\GoogleDrive\HEBEC\Thesis\fig\Spectroscopy';
% f1=stfig('Final plot');
% clf;
% 
% data_colours = [226,75,42;
%                 45,225,245]/255;
% theory_colours = 0.8*data_colours;
% 
% f_pred = opts.const.f_table.g_2_3P_2.(sprintf('e_%s',fmt_name))/1e6;
% fmt_name = strrep(e_level,'^','_');
% for cidx =1:ncat
%     subplot(2,1,cidx)
%     cdata = data.cat{cidx};
%     maxpt = max(cdata.spec.signal);
%     normalizer = 1;
%     B = opts.Bfield(cidx);
%     %        fnum = 7777+cidx;
%            
%     [~,p_cen] = max(cdata.pfits.lor_prms(:,3));
% %     f_offset = cdata.pfits.lor_prms(p_cen,1);
%     f_offset = f_pred;
%     X_data=cdata.spec.freq;
%     all_fit=zeros(1,100);
%     
%     plot(X_data-f_offset,cdata.spec.signal/normalizer,'.','MarkerSize',15,'Color',data_colours(cidx,:))
%     hold on
%     for pidx=1:num_pks
%         if num_pks == 1
%             e_state = opts.e_state;
%         else
%             e_state = opts.e_state{pidx};
%         end
%         e_level = e_state(1:6);
%         fmt_name = strrep(e_level,'^','_');
% %         f_pred = opts.const.f_table.g_2_3P_2.(sprintf('e_%s',fmt_name))/1e6; %MHz
%         f_shift = cdata.zeeman.shift(pidx);
%         cli_header(1,'\nRESULTS for %s transition: field %u\n',e_state,cidx);
%         fprintf('Predicted vacuum freq %.3f MHz\n',f_pred)
% 
%         peak_df = cdata.zeeman.shift(pidx)-cdata.zeeman.shift(p_cen);
%         f_cen = cdata.zeeman.corrected(p_cen);
% 
%         l_mdl = cdata.pfits.lorz{pidx};
%         df = cdata.pfits.lor_prms(pidx,1) - cdata.pfits.lor_prms(p_cen,1);
%         X=cdata.spec.freq;
%         Xfit = linspace(min(X),max(X),100);
%         Yfit = lfun(l_mdl.Coefficients.Estimate,Xfit-cdata.pfits.lor_prms(pidx,1));
%         all_fit = all_fit+Yfit;
%         [ysamp,yci]=predict(l_mdl,(Xfit-cdata.pfits.lor_prms(pidx,1))','Prediction','observation'); 
% 
%         freq_stat_err = cdata.zeeman.shift_unc(pidx) + cdata.pfits.lor_err(pidx,1);
%         % Z shift unc + fit err 
%         theory_error = 0.7+cdata.zeeman.shift_unc(pidx);
%         % Drake's quoted error + zeeman shift error for plot
% 
%         fprintf('Peak %u fitted value:   %.2f(%.2f) MHz\n',pidx,cdata.zeeman.corrected(pidx),freq_stat_err)
%         fprintf('         theory diff:   %.2f MHz\n',cdata.zeeman.corrected(pidx)-f_pred)
%     end
%     if length(opts.e_state)==6
%         X_thry = [-theory_error,-theory_error,theory_error,theory_error]+f_pred+f_shift-f_offset;
%         Y_thry = [0,1.2,1.2,0];
%         thrybox=fill(X_thry,Y_thry,theory_colours(cidx,:),'FaceAlpha',0.3,'EdgeColor',[1,1,1]);
%         uistack(thrybox,'bottom')
%     else
%         f_wrt = 744396208.36;
%         d_stage = zeros(2,5);
%         d_stage(1,:) = [247.8678
%                   -16.2909
%                   -47.3705
%                    16.0309
%                    -9.8264];
%         d_stage(2,:) = [262.8501
%                   -10.1469
%                   -30.3542
%                    11.3518
%                    -8.3259];
%         for thidx = 1:4
%             f_pred = f_wrt;
%             f_shift = d_stage(cidx,thidx+1);
%             X_thry = [-theory_error,-theory_error,theory_error,theory_error]+f_shift;
%             Y_thry = [0,1.2,1.2,0];
%             thrybox=fill(X_thry,Y_thry,theory_colours(cidx,:),'FaceAlpha',0.3,'EdgeAlpha',0);
%             uistack(thrybox,'bottom')
%         end
%     end
%     plot(Xfit-f_offset,all_fit/normalizer,'k:','LineWidth',2)
%     set(gca,'box','off','color','none')
% end
% % set box property to off and remove background color
% 
% % 
% % for thidx = 1:4
% %     X_thry = [-theory_error,-theory_error,theory_error,theory_error]+f_pred+f_shift-f_offset;
% %     Y_thry = [0,1.2,1.2,0];
% %     thrybox=fill(X_thry,Y_thry,theory_colours(cidx,:),'FaceAlpha',0.3,'EdgeColor',[1,1,1]);
% %     uistack(thrybox,'bottom')
% % end
% 
% 
% % wrt 744396208.36
% % stage 1
% %   247.8678
% %   -16.2909
% %   -47.3705
% %    16.0309
% %    -9.8264
% % stage 2
% %   262.8501
% %   -10.1469
% %   -30.3542
% %    11.3518
% %    -8.3259
% 
% 
% % suptitle(sprintf('$%s$',e_level))
% % set(gca,'FontSize',20)
% ylabel('Normalized  loss','FontSize',30)
% % ylim([-0,1.2])
% xlabel(sprintf('$\\nu$-%.2f (MHz)',f_offset),'FontSize',30)
% % suptitle(sprintf('$%s$',opts.e_state))
% % daspect([30,2,1])
% imname = sprintf('out_%s_plot_pretty',e_level);
% filename1 = fullfile(savedir,imname);
%             saveas(f1,[filename1,'.fig']);
%             saveas(f1,[filename1,'.png']);
%             saveas(f1,[filename1,'.eps']);
%             saveas(f1,[filename1,'.svg']);
% cli_header(1,'Done.');
% 
% %% Save to output
% cli_header(0,'Saving output...');
% out_data.data = data.cat;
% out_data.options = opts;
% filename = fullfile(opts.out_dir,'output_and_options.mat');
% save(filename,'out_data','-v7.3')
% fwtext('All Done!')
