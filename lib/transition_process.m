function tr = transition_process(data,opts_tr)

    num_files = numel(data.ai.data);
    ai_time = data.ai.timestamp;
    lv_time = data.lv.time;
    wm_time = data.wm.blue_freq.posix_time;
    
    
    
    %% FETCHING CALIBRATION
    % Can probably combine with preceding loop...
    % calibration changes sign UP when switching to calibration
    % Changes sign DOWN when changing to measurement
    mode_swap = diff(data.lv.calibration);
    % Should guarantee correct partitioning
    first_cals = find(mode_swap > 0);
    last_cals = find(mode_swap < 0);
    if max(last_cals) > num_files
        last_cals = last_cals(last_cals<num_files);
        first_cals = first_cals(first_cals<num_files);
    end
    if max(last_cals) < max(first_cals)
        % Cal switches on, but the return not caught
        last_cals = [last_cals,num_files];
    end
    if min(first_cals) > min(last_cals)
        % The first shots were calibrations but the switch not caught
        first_cals = [1,first_cals];
        % all other cases should be fine?
    end
    num_cal_seg = length(first_cals);
    if length(last_cals) ~= length(first_cals)
        warning('Incomplete calibration segment!')
    end
    if num_cal_seg == 0
        warning('No calibration shots found!')
    end
    
    
    cal_vals = [];
    num_seg_init = max(num_cal_seg,1);
    cal_vals.mean = zeros(num_seg_init,1);
    cal_vals.std = zeros(num_seg_init,1);
    cal_vals.num_shots = zeros(num_seg_init,1);
    cal_vals.time = zeros(num_seg_init,1);
    for ii=1:num_cal_seg
           this_pd_data = data.ai.pd_range(first_cals(ii):last_cals(ii));
           this_time = ai_time(first_cals(ii):last_cals(ii));
           this_pd_data = this_pd_data(~isnan(this_pd_data));
           cal_vals.mean(ii) = nanmean(this_pd_data);
           cal_vals.std(ii)=nanstd(this_pd_data);
           cal_vals.num_shotsum_shots(ii)=length(this_pd_data);           
           cal_vals.time(ii)=nanmean(this_time);
    end
    
    %% Matching timestamps
    sync_shots = [];
    sync_shots.pd_range = zeros(num_files,1);
    sync_shots.wm_set = zeros(num_files,1);
    sync_shots.lv_cal = zeros(num_files,1);
    sync_shots.lv_set = zeros(num_files,1);
    sync_shots.ai_time = zeros(num_files,1);
    sync_shots.wm_time = zeros(num_files,1);
    sync_shots.lv_time = zeros(num_files,1);
    sync_shots.cal_mean = zeros(num_files,1);
    sync_shots.cal_std = zeros(num_files,1);
    sync_shots.cal_N = zeros(num_files,1);
    for ii = 1:num_files
       % Fun job: Write fn/ mod closest_value so it can be vectorized
        this_ai_time = ai_time(ii);
        [this_wm_time,this_wm_idx] = closest_value(wm_time,this_ai_time);
        this_wm_set = data.wm.blue_freq.value(this_wm_idx);
        [this_lv_time,this_lv_idx] = closest_value(lv_time,this_ai_time);
        this_lv_cal = data.lv.calibration(this_lv_idx);
        this_lv_set = data.lv.setpoint(this_lv_idx);
        [~ ,this_cal_idx] = closest_value(cal_vals.time,this_ai_time);
        sync_shots.pd_range(ii) = data.ai.pd_range(ii);
        sync_shots.wm_set(ii) = this_wm_set;
        sync_shots.lv_cal(ii) = this_lv_cal;
        sync_shots.lv_set(ii) = this_lv_set;
        sync_shots.ai_time(ii) = this_ai_time;
        sync_shots.wm_time(ii) = this_wm_time;
        sync_shots.lv_time(ii) = this_lv_time;
        sync_shots.cal_mean(ii) = cal_vals.mean(this_cal_idx);
        sync_shots.cal_std(ii) = cal_vals.std(this_cal_idx);
        sync_shots.cal_N(ii) = cal_vals.num_shots(this_cal_idx);
    end
    
    %Separating test & cal shots
%     ctrl_mask = data.lv.calibration;
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
    sync_msr.calibrated = sync_msr.pd_range-sync_msr.cal_mean; 

    
    % These two lines are inefficient & could be replaced by a single imported item in the
    % interface, but that's a change for later
    all_setpts = unique(data.lv.setpoint);
    all_setpts = all_setpts(~isnan(all_setpts));
    sort_setpt = sort(all_setpts);
    min_set = 2*min(all_setpts);
    max_set = 2*max(all_setpts);

    
    %% Grouping by wavelength
    num_freq_bins = length(all_setpts);
    freq_gap = mean(diff(sort_setpt));
    fbin_edges = linspace(min_set-freq_gap,max_set+freq_gap,num_freq_bins+1)/1e6; %MHz
    [msr_setpts ,msr_set_id]= sort(sync_msr.wm_set);
    freq_stats.freq = zeros(num_freq_bins,1);
    freq_stats.pd_sig = zeros(num_freq_bins,1);
    freq_stats.std_freq = zeros(num_freq_bins,1);
    freq_stats.std_pd = zeros(num_freq_bins,1);
    freq_stats.num_shots = zeros(num_freq_bins,1);
    freq_stats.pd_sig_cal = zeros(num_freq_bins,1);
    freq_stats.pd_err_cal = zeros(num_freq_bins,1);
    for ii=1:num_freq_bins
       fmask_temp = msr_setpts > fbin_edges(ii) & msr_setpts < fbin_edges(ii+1);
       iis = msr_set_id(fmask_temp);
       shots_temp = struct_mask(sync_msr,iis);
       freq_stats.freq(ii) = nanmean(shots_temp.wm_set);
       freq_stats.pd_sig(ii) = nanmean(shots_temp.pd_range);
       freq_stats.freq_err(ii) = nanstd(shots_temp.wm_set);
       freq_stats.pd_sig_err(ii) = nanstd(shots_temp.pd_range);
       freq_stats.num_shots(ii) = sum(fmask_temp);
       freq_stats.pd_sig_cal(ii) = nanmean(shots_temp.pd_range - shots_temp.cal_mean);
    end
    
    % Write output
    tr.shots = sync_shots;
    tr.stats = freq_stats;
    tr.calib = cal_vals;

    
if opts_tr.plot
    %% Compute things for plots


    [cal_hist,cal_bin_edges]= histcounts(sync_cal.pd_range,opts_tr.num_cal_bins);
    cal_bin_cents = 0.5*(cal_bin_edges(1:end-1)+cal_bin_edges(2:end));

    [wm_hist,cal_wm_bin_edges] = histcounts(wm_set_err(wm_msr_mask),opts_tr.num_cal_bins);
    wm_bin_cents = 0.5*(cal_wm_bin_edges(1:end-1)+cal_wm_bin_edges(2:end));
    
    
    mid_setpt = opts_tr.pred_freq*1e3; %MHz
    pd_val_raw = sync_msr.pd_range;
    pd_val_cal = sync_msr.pd_range - sync_msr.cal_mean;
    
    freq_raw_X = sync_msr.wm_set; %MHz
    
    bin_freq_X = freq_stats.freq-mid_setpt; % MHz
    bin_freq_Y = freq_stats.pd_sig_cal;
    bin_freq_Y_err = freq_stats.pd_sig_err./freq_stats.num_shots';
    
    %% Plotting
    sfigure(400);
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


    sfigure(600);
    clf;
    subplot(2,2,1)
    plot((freq_raw_X-mid_setpt),pd_val_raw,'.')
    xlabel('Blue WM reading [MHz]')
    ylabel('PD range')
    title('Raw responses')
    
    
    subplot(2,2,2)
    plot((freq_raw_X-mid_setpt),pd_val_cal,'.')
    xlabel('Blue WM reading [MHz]')
    ylabel('PD range')
    title('Calibrated responses')
    
    subplot(2,2,[3 4])
    errorbar(bin_freq_X,bin_freq_Y,bin_freq_Y_err,'.')
%     xlim([min_set/1e6,max_set/1e6])
    xlabel(sprintf('f-%3.5f [MHz]',mid_setpt))
    ylabel('PD range')
    title('Binned & calibrated response')
    
    suptitle('Processing diagnostics')
    
end
    
    
end