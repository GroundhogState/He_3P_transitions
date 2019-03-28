function tr = bec_transition_process(data,opts_tr)

%     num_files = numel(data.tdc.shot_num);
    tdc_time = data.tdc.time_create_write(:,2);
    ai_time = data.ai.timestamp;
    lv_time = data.lv.time;
    wm_time = data.wm.blue_freq.posix_time;
    
    mid_setpt = opts_tr.pred_freq*1e3; %MHz
    
    start_time = min(lv_time);
    num_files = length(lv_time);
    
    
    %% Matching timestamps
    sync_shots = [];
    sync_shots.num_atoms = zeros(num_files,1);
    sync_shots.wm_set = zeros(num_files,1);
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
       % Fun job: Write fn/ mod closest_value so it can be vectorized
        this_time = lv_time(this_lv_idx); %The LV log is the first thing written, so timestamp to that
        sync_shots.lv_cal(this_lv_idx) = data.lv.calibration(this_lv_idx);
        sync_shots.lv_set(this_lv_idx) = data.lv.setpoint(this_lv_idx);
        sync_shots.lv_time(this_lv_idx) = this_time;
        
        [this_wm_time,this_wm_idx] = closest_value(wm_time,this_time);
        sync_shots.wm_time(this_lv_idx) = this_wm_time;
        sync_shots.wm_set(this_lv_idx) = data.wm.blue_freq.value(this_wm_idx);
        
        [~,this_tdc_idx] = closest_value(tdc_time,this_time);
        this_tdc_idx = this_tdc_idx + 1;
        sync_shots.tdc_time(this_lv_idx) = tdc_time(this_tdc_idx);
        sync_shots.N_atoms(this_lv_idx) = data.tdc.N_atoms(this_tdc_idx)';
        
    end
    
    %% Separating test & cal shots
    ctrl_mask = logical(sync_shots.lv_cal)';
    % Masking for laser setpt errors
    wm_set_err = sync_shots.wm_set - 2*sync_shots.lv_set/1e6; % Error in MHz
    wm_set_mask = ~isnan(sync_shots.wm_set)';
    wm_msr_mask = (abs(wm_set_err) < opts_tr.wm_tolerance)';
    
    
    %Create the masks
    cal_mask = ctrl_mask;
    msr_mask = ~ctrl_mask & wm_set_mask;
    % Break out calibration % measurement blocks
    sync_cal = struct_mask(sync_shots,cal_mask);
    sync_msr = struct_mask(sync_shots,msr_mask);
    
   
    %% create a calibration model
    % v1: Interpolation. Fast, rough.
    interp_mdl = interp1(sync_cal.tdc_time,sync_cal.N_atoms,sync_msr.tdc_time);
    sync_msr.interp_mdl = interp_mdl;
       

    % These lines are inefficient & could be replaced by a single imported item in the
    % interface, but that's a change for later
    all_setpts = unique(data.lv.setpoint);
    all_setpts = all_setpts(~isnan(all_setpts));
    sort_setpt = sort(all_setpts);
    min_set = 2*min(all_setpts);
    max_set = 2*max(all_setpts);

    
    %% Grouping by wavelength (Independent variable)
    num_freq_bins = length(all_setpts);
    freq_gap = mean(diff(sort_setpt));
    fbin_edges = linspace(min_set-freq_gap,max_set+freq_gap,num_freq_bins+1)/1e6; %MHz
    [msr_setpts ,msr_set_id]= sort(sync_msr.wm_set);
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
           freq_stats.freq(ii) = nanmean(shots_temp.wm_set);
           freq_stats.sig(ii) = nanmean(shots_temp.N_atoms);
           freq_stats.freq_err(ii) = nanstd(shots_temp.wm_set);
           freq_stats.num_shots(ii) = sum(fmask_temp);
           freq_stats.sig_cal(ii) = nanmean(shots_temp.N_atoms./shots_temp.interp_mdl');
           freq_stats.sig_err(ii) = nanstd(shots_temp.N_atoms./shots_temp.interp_mdl');
           freq_stat_mask(ii) = abs(freq_stats.freq(ii)) > 0;
       end
    end
    freq_stats = struct_mask(freq_stats,logical(freq_stat_mask));
    
    
    %% Fitting the data
%     inv_peak = 1-freq_stats.sig_cal;
    
%     [peak_guess,peak_idx]= min(freq_stats.sig_cal);
%     var_guess = 15;
%     % Fit a dumb gaussian
%     mdlfun = @(p,x) 1 - p(1)*exp(-0.5*((x-p(2))/p(3)).^2);
%     beta0 = [1-peak_guess,freq_stats.freq(peak_idx),var_guess];
%     fit = fitnlm(freq_stats.freq,freq_stats.sig_cal,mdlfun,beta0);
    
    
    
    
    %% Write output
    tr.sense = sync_msr;
    tr.calib = sync_cal;
    tr.stats = freq_stats;


if opts_tr.plot
    %% Plotting
    % % Compute things

    [cal_hist,cal_bin_edges]= histcounts(sync_cal.N_atoms,opts_tr.num_cal_bins);
    cal_bin_cents = 0.5*(cal_bin_edges(1:end-1)+cal_bin_edges(2:end));

    [wm_hist,cal_wm_bin_edges] = histcounts(wm_set_err(wm_msr_mask),opts_tr.num_cal_bins);
    wm_bin_cents = 0.5*(cal_wm_bin_edges(1:end-1)+cal_wm_bin_edges(2:end));
    

  
    plot_raw_T = sync_msr.tdc_time-start_time;
    plot_raw_X = sync_msr.wm_set-mid_setpt;
    plot_raw_Y = sync_msr.N_atoms;
    
    plot_cal_T = sync_cal.tdc_time-start_time;
    plot_cal_X = sync_cal.wm_set-mid_setpt;
    plot_cal_Y = sync_cal.N_atoms;
    
    plot_mdl_Y = sync_msr.N_atoms./interp_mdl';
    
    plot_sig_X = freq_stats.freq-mid_setpt;
    plot_sig_Y = freq_stats.sig_cal;
    plot_sig_Y_err = freq_stats.sig_err;
    plot_sig_Y_err = plot_sig_Y_err.*~isnan(plot_sig_Y_err);
    
    
    
    %% Plotting
    % Plotting the calibration model
        
    f1=sfigure(5000);
    clf;
    subplot(2,1,1)
    plot(plot_raw_T,plot_raw_Y,'.')
    hold on
    plot(plot_cal_T,plot_cal_Y,'.')
    plot(plot_raw_T,interp_mdl)
    xlabel('Time elapsed')
    ylabel('Counts')
    legend('Measurement shots','calibration shots','Model')
    title('Raw data')

    subplot(2,1,2)
    plot(plot_raw_X,plot_mdl_Y,'.')
    xlabel(sprintf('f-%3.5f [MHz]',mid_setpt))
    ylabel('N ratio')
    title('Model-calibrated signal')    
    suptitle('Calibration model')
    
    % Plotting the histograms
    f2=sfigure(500);
    clf
    subplot(2,2,1)
    area(cal_bin_cents,cal_hist)
    xlabel('Value [V]')
    ylabel('Counts')
    title(sprintf('Range calibration, N=%u',sum(cal_mask)))

    subplot(2,2,2)
    plot(sync_shots.wm_set,'x')
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
    plot(ai_time,'.')
    xlabel('Shot number')
    ylabel('Time')
    title('Recorded timestamps')
    legend('WM','TDC','LV','AI')

    % Plotting the result
    f3=sfigure(501);
    clf;
    errorbar(plot_sig_X,plot_sig_Y,plot_sig_Y_err,'.')
    hold on
%     plot(freq_stats.freq-mid_setpt,mdlfun(beta0,freq_stats.freq),'kx')
%     plot(freq_stats.freq-mid_setpt,mdlfun(fit.Coefficients.Estimate,freq_stats.freq),'ro')
%     plot(plot_sig_X,mdlfun(fit.Coefficients.Estimate,plot_sig_X+mid_setpt))
    xlabel(sprintf('f-%3.5f [MHz]',mid_setpt))
    ylabel('Atom Number Ratio')

    suptitle('Binned absorption spectrum')
    
    filename1 = fullfile(opts_tr.out_dir,sprintf('%s_raw_calib',mfilename));
    saveas(f1,[filename1,'.fig']);
    saveas(f1,[filename1,'.png'])

    filename2 = fullfile(opts_tr.out_dir,sprintf('%s_diagnostic',mfilename));
    saveas(f2,[filename2,'.fig']);
    saveas(f2,[filename2,'.png'])

    filename3 = fullfile(opts_tr.out_dir,sprintf('%s_spectrum',mfilename));
    saveas(f3,[filename3,'.fig']);
    saveas(f3,[filename3,'.png'])
    
end

end