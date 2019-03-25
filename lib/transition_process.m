function tr = transition_process(data,opts_tr)


    

%     if ~isnan(opts_tr.num_files)
    num_files = numel(data.ai.file_names);

    % MATCHING TIMESTAMPS
    ai_time = data.ai.timestamp;
    lv_time = data.lv.time;
    % Needs moar logging?!
    wm_time = data.wm.blue_freq.posix_time;
    
    
    
    % % FETCHING CALIBRATION
    % Can probably combine with preceding loop...
    % calibration changes sign UP when switching to calibration
    % Changes sign DOWN when changing to measurement
    mode_swap = diff(data.lv.calibration);
    % Should guarantee correct partitioning
    first_cals = find(mode_swap > 0);
    last_cals = find(mode_swap < 0);
    if max(last_cals) < max(first_cals)
        % Cal switches on, but the return not caught
        last_cals = [last_cals,num_files];
    end
    if min(first_cals) > min(last_cals)
        % The first shots were calibrations but the switch not caught
        first_cals = [1,first_cals];
        % all other cases should be fine?
    end
    if length(last_cals) ~= length(first_cals)
        warning('Incomplete calibration segment!')
    end
    num_cal_seg = length(last_cals);
    cal_vals = zeros(num_cal_seg,4);
    for ii=1:num_cal_seg
       this_data = data.ai.pd_range(first_cals(ii):last_cals(ii));
       this_time = ai_time(first_cals(ii):last_cals(ii));
       this_data = this_data(~isnan(this_data));
       cal_vals(ii,:) = [nanmean(this_data),nanstd(this_data),length(this_data),nanmean(this_time)];
    end
    
    
    sync_shots = zeros(num_files,10);
%     cal_mask = zeros(num_files,1);
    for shot_idx = 1:num_files
       % Fun job: Write fn/ mod closest_value so it can be vectorized
        this_ai_time = ai_time(shot_idx);
        [this_wm_time,this_wm_idx] = closest_value(wm_time,this_ai_time);
        this_wm_set = data.wm.blue_freq.value(this_wm_idx);
        [this_lv_time,this_lv_idx] = closest_value(lv_time,this_ai_time);
        this_lv_cal = data.lv.calibration(this_lv_idx);
        this_lv_set = data.lv.setpoint(this_lv_idx);
        [this_cal_time ,this_cal_time_idx] = closest_value(cal_vals(:,4),this_ai_time);
        this_cal = cal_vals(this_cal_time_idx,1:3);
%         cal_mask(shot_idx) = this_lv_cal;
        %                           pd_range                %act freq      %cal     % setpt
        this_sync = [data.ai.pd_range(shot_idx),this_wm_set,this_lv_cal,this_lv_set,ai_time(shot_idx),this_wm_time,this_lv_time];
        sync_shots(shot_idx,:) = [this_sync,this_cal];
    end
    

    
    %Separating test & cal shots
%     ctrl_mask = data.lv.calibration;
    ctrl_mask = sync_shots(:,3)';
    % Masking for laser setpt errors
    wm_set_err = sync_shots(:,2) - 2*sync_shots(:,4)/1e6; % Error in MHz
    wm_set_mask = ~isnan(sync_shots(:,2))';
    wm_msr_mask = (abs(wm_set_err) < opts_tr.wm_tolerance)';
    
    %ctrl_mask = logical([ctrl_mask,0,0]);
    
    %Create the masks
    cal_mask = logical(ctrl_mask);
    msr_mask = ~ctrl_mask & wm_set_mask;
    
    sync_cal = sync_shots(cal_mask,:);
    sync_msr = sync_shots(msr_mask,:);
    sync_msr = [sync_msr,sync_msr(:,1)-sync_msr(:,8)]; 

    
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
    [msr_setpts ,msr_set_id]= sort(sync_msr(:,2));
    freq_stats = zeros(num_freq_bins,7);
    for ii=1:num_freq_bins
%        fmask_temp = fast_sorted_mask(msr_setpts,fbin_edges(ii),fbin_edges(ii+1));
       fmask_temp = msr_setpts > fbin_edges(ii) & msr_setpts < fbin_edges(ii+1);
       shot_idxs = msr_set_id(fmask_temp);
       shots_temp = sync_msr(shot_idxs,:);
       % mean freq,     mean sig,       freq std,       sig std         num shots
       freq_stats(ii,:) = [nanmean(shots_temp(:,2)),nanmean(shots_temp(:,1)),nanstd(shots_temp(:,2)),nanstd(shots_temp(:,1)),sum(fmask_temp),nanmean(shots_temp(:,11)),nanstd(shots_temp(:,11))];
    end
    
    % Write output
    tr.shots = sync_shots;
    tr.stats = freq_stats;
    tr.calib = cal_vals;

    
    
if opts_tr.plot
    %% Compute things for plots


    [cal_hist,cal_bin_edges]= histcounts(sync_cal(:,1),opts_tr.num_cal_bins);
    cal_bin_cents = 0.5*(cal_bin_edges(1:end-1)+cal_bin_edges(2:end));

    [wm_hist,cal_wm_bin_edges] = histcounts(wm_set_err(wm_msr_mask),opts_tr.num_cal_bins);
    wm_bin_cents = 0.5*(cal_wm_bin_edges(1:end-1)+cal_wm_bin_edges(2:end));
    
    
    num_shots_final = 
    
    mid_setpt = nanmedian(sync_msr(:,2));
    pd_val_raw = sync_msr(:,1);
    pd_val_cal = sync_msr(:,1) - sync_msr(:,8);
    
    bin_freq_X = freq_stats(:,1)-mid_setpt;
    bin_freq_Y = freq_stats(:,6);
    bin_freq_Y_err = freq_stats(:,7)./freq_stats(:,5);
    
    %% Plotting
    sfigure(400);
    clf


    subplot(2,2,1)
    area(cal_bin_cents,cal_hist)
    xlabel('Value [V]')
    ylabel('Counts')
    title(sprintf('Range calibration, N=%u',sum(cal_mask)))

    subplot(2,2,2)
    plot(sync_shots(:,2),'x')
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
    plot((sync_msr(:,2)-mid_setpt),pd_val_raw,'.')
    xlabel('Blue WM reading [MHz]')
    ylabel('PD range')
    title('Raw responses')
    
    
    subplot(2,2,2)
    plot((sync_msr(:,2)-mid_setpt),pd_val_cal,'.')
    xlabel('Blue WM reading [MHz]')
    ylabel('PD range')
    title('Calibrated responses')
    
    subplot(2,2,[3 4])
    errorbar(bin_freq_X,bin_freq_Y,bin_freq_Y_err,'.')
%     xlim([min_set/1e6,max_set/1e6])
    xlabel(sprintf('f-%3.5f [MHz]',mid_setpt))
    ylabel('PD range')
    title(sprintf('Binned & calibrated responses, %u shots',num_shots_final))
    
    suptitle('Processing diagnostics')
    
    
    sfigure(601)
    plot(bin_freq_X,log(abs(bin_freq_Y)))
%     xlim([min_set/1e6,max_set/1e6])
    xlabel(sprintf('f-%3.5f [MHz]',mid_setpt))
    ylabel('PD range')
    title('Binned & calibrated response')
    
end
    
    
end