function tr = bec_transition_process(data,opts_tr)


    %% %%%%%%%%%%%% PRECONDITIONING THE DATA 
%     num_files = numel(data.tdc.shot_num);
    tdc_time = data.tdc.time_create_write(:,2);
%     ai_time = data.ai.timestamp;
    lv_time = data.lv.time;
    wm_time = data.wm.feedback.posix_time;
    
    all_setpts = unique(data.lv.setpoint);
    lv_set_range = [min(all_setpts)/1e6,max(all_setpts)/1e6]; %MHz
    mid_setpt = mean(lv_set_range); %MHz
    
    start_time = min(lv_time);
    num_files = length(lv_time);
    
    
    %% Matching timestamps
    % Should:   Truncate all logs by time to fit within smallest logged
    % window
    % Plot all diagnostics
    % Extract as function/consider generalizing
    % What's with the 10MHz offset?!
    sync_shots = [];
    sync_shots.num_atoms = zeros(num_files,1);
    sync_shots.wm_setpt = zeros(num_files,1);
    sync_shots.lv_cal = zeros(num_files,1);
    sync_shots.lv_set = zeros(num_files,1);
    sync_shots.tdc_time = zeros(num_files,1);
    sync_shots.wm_time = zeros(num_files,1);
    sync_shots.lv_time = zeros(num_files,1);
    sync_shots.cal_mean = zeros(num_files,1);
    sync_shots.cal_std = zeros(num_files,1);
    sync_shots.cal_N = zeros(num_files,1);
    sync_shots.cal_interp = zeros(num_files,1);
    for this_lv_idx = 1:num_files
       % Fun job: Write fn/ modify closest_value so it can be vectorized
        this_time = lv_time(this_lv_idx); %The LV log is the first thing written, so timestamp to that
        sync_shots.lv_cal(this_lv_idx) = data.lv.calibration(this_lv_idx);
        sync_shots.lv_set(this_lv_idx) = data.lv.setpoint(this_lv_idx);
        sync_shots.lv_time(this_lv_idx) = this_time;
        
        [this_wm_time,this_wm_idx] = closest_value(wm_time,this_time+2);
        sync_shots.wm_time(this_lv_idx) = this_wm_time;
        sync_shots.wm_setpt(this_lv_idx) = data.wm.feedback.actual(this_wm_idx);
%         sync_shots.wm_blue = data.wm.blue_freq(this_wm_idx);
        
        % CORRECT FOR AOM OFFSET
        sync_shots.probe_set(this_lv_idx) = 2*sync_shots.wm_setpt(this_lv_idx) - opts_tr.aom_freq;
        
        [~,this_tdc_idx] = closest_value(tdc_time,this_time);
        this_tdc_idx = this_tdc_idx + 1;
        sync_shots.tdc_time(this_lv_idx) = tdc_time(this_tdc_idx);
        sync_shots.N_atoms(this_lv_idx) = data.tdc.N_atoms(this_tdc_idx)';
    end
    


    
    %% Separating test & cal shots
    ctrl_mask = logical(sync_shots.lv_cal)';
    % Masking for laser setpt errors
    wm_set_err = sync_shots.wm_setpt - 2*sync_shots.lv_set/1e6; % Error in MHz
    wm_set_mask = ~isnan(sync_shots.wm_setpt)';
    wm_msr_mask = (abs(wm_set_err) < opts_tr.wm_tolerance)';
    
    
    %Create the masks
    cal_mask = ctrl_mask;
    msr_mask = ~ctrl_mask & wm_set_mask;
    % Break out calibration % measurement blocks
    sync_cal = struct_mask(sync_shots,cal_mask);
    sync_msr = struct_mask(sync_shots,msr_mask);
    
    lv_mask = ~isnan(data.lv.setpoint);
    sfigure(502);
    clf;
    subplot(2,2,1)
    plot(sync_shots.lv_time-start_time,sync_shots.lv_set/1e6,'x');
    hold on
    plot(wm_time-start_time,data.wm.feedback.actual,'x');
    legend('2*LV set point','WM blue setpt')  
    xlabel('Time elapsed')
    title('Raw channel inputs')
    subplot(2,2,2)
    plot(sync_shots.lv_time-start_time,sync_shots.lv_set/1e6,'x');
    hold on
    plot(sync_shots.lv_time-start_time,sync_shots.wm_setpt,'x');
    legend('2*LV set point','WM blue setpt')  
    xlabel('Time elapsed')
    ylabel(sprintf('Frequency - %3f',mid_setpt))
    legend('LabView log set','WM record')
    title('Time-matched setpoints')
    
    subplot(2,2,3)
    plot(sync_msr.lv_time,sync_msr.lv_set/1e6 - sync_msr.wm_setpt,'x')
    xlabel('Time elapsed')
    ylabel('2*LV set - WM set')
    title('Setpoint error')
    
    subplot(2,2,4)
    plot(sync_msr.lv_time,sync_msr.probe_set,'x')
    
    title('AOM-corrected scanned range')
    
    suptitle('Probe beam setpoint matching')
    
   
    %% %%%%%%%%% COMPUTATION ON THE DATA
    
    %% create a calibration model
    % v1: Interpolation. Fast, rough.
    interp_mdl = interp1(sync_cal.tdc_time,sync_cal.N_atoms,sync_msr.tdc_time,'spline');
    sync_msr.interp_mdl = interp_mdl;
       
    min_set = min(sync_msr.probe_set);%2*min(all_setpts);
    max_set = max(sync_msr.probe_set);%2*max(all_setpts);


    
    %% Grouping by wavelength (Independent variable)
    num_freq_bins = round(range(all_setpts/1e6));
    if ~isfield(opts_tr,'num_freq_bins')
        if ~isnan(opts_tr.num_freq_bins) 
            num_freq_bins = opts_tr.num_freq_bins;
        end
    end
    freq_gap = 0;%mean(diff(sort_setpt));
    fbin_edges = linspace(min_set-freq_gap,max_set+freq_gap,num_freq_bins+1);%/1e6 %MHz
    [msr_setpts ,msr_set_id]= sort(sync_msr.probe_set);
    freq_stats.freq = zeros(num_freq_bins,1);
    freq_stats.sig = zeros(num_freq_bins,1);
    freq_stats.std_freq = zeros(num_freq_bins,1);
    freq_stats.std_pd = zeros(num_freq_bins,1);
    freq_stats.num_shots = zeros(num_freq_bins,1);
    freq_stats.sig_cal = zeros(num_freq_bins,1);
    freq_stats.pd_err_cal = zeros(num_freq_bins,1);
    freq_stat_mask = zeros(num_freq_bins,1);
    for ii=1:num_freq_bins
       fmask_temp = msr_setpts > fbin_edges(ii) & msr_setpts < fbin_edges(ii+1);
       if sum(fmask_temp) >0
           iis = msr_set_id(fmask_temp);
           shots_temp = struct_mask(sync_msr,iis);
           freq_stats.freq(ii) = nanmean(shots_temp.probe_set);
           freq_stats.sig(ii) = nanmean(shots_temp.N_atoms);
           freq_stats.freq_err(ii) = nanstd(shots_temp.probe_set);
           freq_stats.num_shots(ii) = sum(fmask_temp);
           freq_stats.sig_cal(ii) = nanmean(shots_temp.N_atoms./shots_temp.interp_mdl');
           freq_stats.sig_std(ii) = nanstd(shots_temp.N_atoms./shots_temp.interp_mdl');
           freq_stats.sig_err(ii) = freq_stats.sig_std(ii)/sqrt(freq_stats.num_shots(ii));
           freq_stat_mask(ii) = abs(freq_stats.freq(ii)) > 0;
       end
    end
    freq_stats = struct_mask(freq_stats,logical(freq_stat_mask));
    
    
    %% Fitting the data
   
    [peak_guess,peak_idx]= min(freq_stats.sig_cal);
    peak_freq = freq_stats.freq(peak_idx);
    freqs_centred_fit = sync_msr.probe_set-peak_freq;
    var_guess = 15;
    % Fit a gaussian
    gaus_mdlfun = @(p,x) 1 - p(1)*exp(-0.5*((x-p(2))/p(3)).^2);
    beta0 = [1-peak_guess,0,var_guess];
    fit_gaus = fitnlm(freqs_centred_fit,sync_msr.N_atoms./interp_mdl',gaus_mdlfun,beta0);

    % Fit Lorentzian
    var_guess = 15;
    lor_mdlfun = @(b,x) b(1)./((x-b(2)).^2+b(3))+1;
    beta0 = [-peak_guess,0,var_guess];
    fit_lor = fitnlm(freqs_centred_fit,sync_msr.N_atoms./interp_mdl',lor_mdlfun,beta0);
    lor_prms = fit_lor.Coefficients.Estimate;
    % Voigt profile?
    
    
    
    %extract useful parameters
    lorz_cen = fit_lor.Coefficients.Estimate(2)+peak_freq;
    lorz_unc = fit_lor.Coefficients.SE(2);
    lorz_width = 2*sqrt(abs(fit_lor.Coefficients.Estimate(3)));
    lorz_width_unc = 2*fit_lor.Coefficients.Estimate(3);
    lorz_wav = (299792458/lorz_cen)*1e3;
    lorz_cen_Hz = lorz_cen*1e6;
    lorz_unc_Hz = lorz_unc*1e6;
    lorz_wav_unc = lorz_unc_Hz*299792458/((lorz_cen_Hz)^2) * 1e9; % convert back to nm
    

    gaus_cen = fit_gaus.Coefficients.Estimate(2)+peak_freq;%+opts_tr.aom_freq;
    gaus_unc = fit_gaus.Coefficients.SE(2);
    gaus_width = fit_gaus.Coefficients.Estimate(3);
    gaus_width_unc = fit_gaus.Coefficients.SE(3);
    gaus_cen_Hz = lorz_cen*1e6;
    gaus_unc_Hz = lorz_unc*1e6;
    gaus_wav = (299792458/gaus_cen)*1e3;
    gaus_wav_unc = gaus_unc_Hz*299792458/((gaus_cen_Hz)^2) * 1e9; % convert back to nm
    
    % 
    freqs_centred = sync_msr.probe_set-lorz_cen;
    
    gaus_prms = fit_gaus.Coefficients.Estimate;
    gaus_cen_prms = gaus_prms;
    gaus_cen_prms(2) = gaus_cen_prms(2) - fit_lor.Coefficients.Estimate(2);
    lor_cen_prms = lor_prms;
    lor_cen_prms(2) = 0;
    
    
    
    %% Write output
    tr.sense = sync_msr;
    tr.calib = sync_cal;
    tr.stats = freq_stats;
    tr.gaus_fit = [gaus_cen,gaus_unc,gaus_width,gaus_width_unc,gaus_wav,gaus_wav_unc];
    tr.lorz_fit = [lorz_cen,lorz_unc,lorz_width,lorz_width_unc,lorz_wav,lorz_wav_unc];



if opts_tr.plot
    %% Plotting
    % % Compute things

    [cal_hist,cal_bin_edges]= histcounts(sync_cal.N_atoms,opts_tr.num_cal_bins);
    cal_bin_cents = 0.5*(cal_bin_edges(1:end-1)+cal_bin_edges(2:end));

    [wm_hist,cal_wm_bin_edges] = histcounts(wm_set_err(wm_msr_mask),opts_tr.num_cal_bins);
    wm_bin_cents = 0.5*(cal_wm_bin_edges(1:end-1)+cal_wm_bin_edges(2:end));
  
%     plot_cen = lorz_cen - peak_freq;
    plot_pred = opts_tr.pred_freq.*1e3-lorz_cen;
%     plot_pred_range = 
    
    plot_set_X = 2*sort(all_setpts);
    plot_set_X = plot_set_X(~isnan(plot_set_X))/1e6; %MHz
    
    plot_raw_T = sync_msr.tdc_time-start_time;
    plot_raw_X = freqs_centred;
    plot_raw_Y = sync_msr.N_atoms;
    
    plot_fit_X = linspace(min([freqs_centred,plot_pred-1]),max([freqs_centred,plot_pred+1]),1e4);
    plot_fit_Y_lor = lor_mdlfun(lor_cen_prms,plot_fit_X);
    plot_fit_Y_gaus = gaus_mdlfun(gaus_cen_prms,plot_fit_X);
    
    plot_cal_T = sync_cal.tdc_time-start_time;
    plot_cal_X = sync_cal.probe_set;
    plot_cal_Y = sync_cal.N_atoms;
    
    plot_mdl_Y = sync_msr.N_atoms./interp_mdl';
    
    plot_sig_X = freq_stats.freq-lorz_cen;
    plot_sig_Y = freq_stats.sig_cal;
    plot_sig_Y_err = freq_stats.sig_err;
    plot_sig_Y_err = plot_sig_Y_err.*~isnan(plot_sig_Y_err);
    
    
    
    %% Plotting
    % Plotting the diagnostics
    f2=sfigure(500);
    clf
    subplot(2,2,1)
    area(cal_bin_cents,cal_hist)
    xlabel('Value [V]')
    ylabel('Counts')
    title(sprintf('Range calibration, N=%u',sum(cal_mask)))

    subplot(2,2,2)
    plot(sync_shots.wm_setpt,'x')
    xlabel('Shot number')
    ylabel('Measured WM freq')
    title('Probe setpt')
    
    subplot(2,2,3)
    area(wm_bin_cents,wm_hist)
    xlabel('Value [MHz]')
    ylabel('Counts')
    title(sprintf('Laser setpoint error, N=%u',length(wm_set_err(wm_set_mask))))
    
    subplot(2,2,4)
    plot(wm_time,'.')
    hold on
    plot(tdc_time,'.')
    plot(lv_time,'.')
    xlabel('Shot number')
    ylabel('Time')
    title('Recorded timestamps')
    legend('WM','TDC','LV')
    

    
    
    % Plotting the calibration model
        
    f1=sfigure(5000);
    clf;
    subplot(4,1,1)
    plot(plot_raw_T,plot_raw_Y,'.')
    hold on
    plot(plot_cal_T,plot_cal_Y,'.')
    plot(plot_raw_T,interp_mdl)
    xlabel('Time elapsed')
    ylabel('Counts')
    legend('Measurement shots','calibration shots','Model')
    title('Raw data')

    subplot(4,1,2)
    plot(plot_raw_T,plot_mdl_Y,'.')
    xlabel('Time elapsed')
    ylabel('N ratio')
    title('Model-calibrated signal')    
    suptitle('Calibration model')
    
    subplot(4,1,3)
    plot(sync_shots.lv_time-start_time,sync_shots.wm_setpt,'.')
    xlabel('Time elapsed')
    ylabel('Raw WM set point')
    
    
    subplot(4,1,4)
    plot(plot_raw_X,plot_mdl_Y,'.')
    hold on
    plot(plot_fit_X,plot_fit_Y_lor,'r-')
    plot(plot_fit_X,plot_fit_Y_gaus,'k--')
    xlabel(sprintf('f-%3.5f [MHz]'))
    ylabel('N ratio')
    title('Model-calibrated signal')    
    suptitle('Calibration model')
    legend('Calibrated data','Lorentzian fit','Gaussian fit')
    


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
    xlabel(sprintf('f - %6f (MHz)',lorz_cen))
    ylabel('Atom Number Ratio')
    legend('Frequency-binned data','Lorentzian fit','Gaussian fit','Theory Value','Location','Best')

    suptitle([opts_tr.tr_name,sprintf(' absorption peak @ %.2f±(%.3f) MHz, FWHM~%.2fMHz',lorz_cen,lorz_unc,lorz_width)])
    title(sprintf('Vacuum wavelength %.10f nm ± %.6f fm, 2*%u shots',lorz_wav,1e6*lorz_wav_unc,sum(msr_mask)))
    
    
    % Saving the figures
    filename1 = fullfile(opts_tr.out_dir,sprintf('%s_raw_calib',mfilename));
    saveas(f1,[filename1,'.fig']);
    saveas(f1,[filename1,'.png'])

    filename2 = fullfile(opts_tr.out_dir,sprintf('%s_diagnostic',mfilename));
    saveas(f2,[filename2,'.fig']);
    saveas(f2,[filename2,'.png'])

    filename3 = fullfile(opts_tr.out_dir,sprintf('%s_spectrum',opts_tr.tr_name));
    saveas(f3,[filename3,'.fig']);
    saveas(f3,[filename3,'.png'])
    
end

end