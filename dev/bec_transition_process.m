function tr = bec_transition_process(data,opts_tr)

    num_files = numel(data.tdc.shot_num);
    tdc_time = data.tdc.time_create_write;
    ai_time = data.ai.timestamp;
    lv_time = data.lv.time;
    wm_time = data.wm.blue_freq.posix_time;
   
    
    %% FETCHING CALIBRATION
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
           this_data = data.tdc.N_atoms(first_cals(ii):last_cals(ii));
           this_time = tdc_time(first_cals(ii):last_cals(ii));
           cal_vals.mean(ii) = nanmean(this_data);
           cal_vals.std(ii)=nanstd(this_data);
           cal_vals.num_shots(ii)=length(this_data);           
           cal_vals.time(ii)=nanmean(this_time);
    end
    %% create a calibration model
%     anal_opts.cal_mdl.smooth_time=100;
%     anal_opts.cal_mdl.plot=true;
%     anal_opts.cal_mdl.global=anal_opts.global;
%     data.cal=make_cal_model(anal_opts.cal_mdl,data);



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
    for ii = 1:num_files
       % Fun job: Write fn/ mod closest_value so it can be vectorized
        this_time = tdc_time(ii);
        [this_wm_time,this_wm_idx] = closest_value(wm_time,this_time);
        this_wm_set = data.wm.blue_freq.value(this_wm_idx);
        [~,this_lv_idx] = closest_value(lv_time,this_time);
%         this_lv_idx = this_lv_idx - 1;
        this_lv_time = data.lv.time(this_lv_idx);
        this_lv_cal = data.lv.calibration(this_lv_idx);
        this_lv_set = data.lv.setpoint(this_lv_idx);
        [~ ,this_cal_idx] = closest_value(cal_vals.time,this_time);
        sync_shots.N_atoms(ii) = data.tdc.N_atoms(ii)';
        sync_shots.wm_set(ii) = this_wm_set;
        sync_shots.lv_cal(ii) = this_lv_cal;
        sync_shots.lv_set(ii) = this_lv_set;
        sync_shots.ai_time(ii) = this_time;
        sync_shots.wm_time(ii) = this_wm_time;
        sync_shots.lv_time(ii) = this_lv_time;
        sync_shots.cal_mean(ii) = cal_vals.mean(this_cal_idx);
        sync_shots.cal_std(ii) = cal_vals.std(this_cal_idx);
        sync_shots.cal_N(ii) = cal_vals.num_shots(this_cal_idx);
    end
    
    %Separating test & cal shots
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
           if num_cal_seg ==0
               shots_temp.cal_mean = 1;
           end
           freq_stats.sig_cal(ii) = nanmean(shots_temp.N_atoms ./ shots_temp.cal_mean');
           freq_stats.sig_err(ii) = nanstd(shots_temp.N_atoms./ shots_temp.cal_mean');
           freq_stat_mask(ii) = abs(freq_stats.freq(ii)) > 0;
       end
    end
    freq_stats = struct_mask(freq_stats,logical(freq_stat_mask));
    
    %% Write output
    tr.shots = sync_shots;
    tr.stats = freq_stats;
    tr.calib = cal_vals;


if opts_tr.plot
    %% Plotting
    % % Compute things

    [cal_hist,cal_bin_edges]= histcounts(sync_cal.N_atoms,opts_tr.num_cal_bins);
    cal_bin_cents = 0.5*(cal_bin_edges(1:end-1)+cal_bin_edges(2:end));

    [wm_hist,cal_wm_bin_edges] = histcounts(wm_set_err(wm_msr_mask),opts_tr.num_cal_bins);
    wm_bin_cents = 0.5*(cal_wm_bin_edges(1:end-1)+cal_wm_bin_edges(2:end));
    
    
    mid_setpt = opts_tr.pred_freq*1e3; %MHz
    plot_raw_X = sync_msr.wm_set-mid_setpt;
    plot_raw_Y = sync_msr.N_atoms;
    plot_cal_X = sync_cal.wm_set-mid_setpt;
    plot_cal_Y = sync_cal.N_atoms;
    plot_sig_X = freq_stats.freq-mid_setpt;
    plot_sig_Y = freq_stats.sig_cal;
    plot_sig_Y_err = freq_stats.sig_err;
    plot_sig_Y_err = plot_sig_Y_err.*~isnan(plot_sig_Y_err);

    
    %% Fitting
       
    
    
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


    sfigure(601);
    clf;
    subplot(2,1,1)
    plot(plot_raw_X,plot_raw_Y,'.')
    hold on
    plot(plot_cal_X,plot_cal_Y,'.')
    xlabel(sprintf('f-%3.5f [MHz]',mid_setpt))
    ylabel('Atom Number')
    legend('Calibration','Interrogation')
    title('Raw')
       
    subplot(2,1,2)
    errorbar(plot_sig_X,plot_sig_Y,plot_sig_Y_err,'.')
    xlabel(sprintf('f-%3.5f [MHz]',mid_setpt))
    ylabel('Atom Number Ratio')
    if num_cal_seg >0
        title('Binned & Calibrated')
    elseif num_cal_seg == 0
        title('Binned, UNCALIBRATED')
    end
    
    suptitle('Processing diagnostics')
    
end

end