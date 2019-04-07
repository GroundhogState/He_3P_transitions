function sync_data = match_timestamps(data,opts)

    %% Matching timestamps
    % Should:   Truncate all logs by time to fit within smallest logged
    % window
    % Plot all diagnostics
    % Extract as function/consider generalizing
    % What's with the 10MHz offset?!
    
    header({0,'Correlating timestamps...'})
    
    tdc_time = data.tdc.time_create_write(:,2);
%     ai_time = data.ai.timestamp;
    lv_time = data.lv.time;
    wm_time = data.wm.feedback.posix_time;
    file_idxs = 1:length(tdc_time);
    
    start_time = min(lv_time);
    all_setpts = unique(data.lv.setpoint);
    lv_set_range = [min(all_setpts)/1e6,max(all_setpts)/1e6]; %MHz
    mid_setpt = mean(lv_set_range); 
    
    if isfield(opts,'ignorefiles')
       warning('Manually ignoring %u of %u files!!! Check opts.tr.filemask',length(opts.ignorefiles),length(file_idxs))
       file_idxs = setdiff(file_idxs,opts.ignorefiles);
       num_files = length(file_idxs); %filemask is a list of shots to drop
    else
        num_files = length(file_idxs);
    end
    
    
    sync_shots = [];
    sync_shots.wm_setpt = zeros(num_files,1);
    sync_shots.lv_cal = zeros(num_files,1);
    sync_shots.lv_set = zeros(num_files,1);
    sync_shots.tdc_time = zeros(num_files,1);
    sync_shots.wm_time = zeros(num_files,1);
    sync_shots.lv_time = zeros(num_files,1);

    for idx = 1:num_files
        this_lv_idx = file_idxs(idx); % To allow for manual masking
       % Fun job: Write fn/ modify closest_value so it can be vectorized
        this_time = lv_time(this_lv_idx); %The LV log is the first thing written, so timestamp to that
        sync_shots.lv_cal(idx) = data.lv.calibration(this_lv_idx);
        sync_shots.lv_set(idx) = data.lv.setpoint(this_lv_idx);
        sync_shots.lv_time(idx) = this_time;
        
        [this_wm_time,this_wm_idx] = closest_value(wm_time,this_time+2);
        sync_shots.wm_time(idx) = this_wm_time;
        sync_shots.wm_setpt(idx) = data.wm.feedback.actual(this_wm_idx);
%         sync_shots.wm_blue = data.wm.blue_freq(this_wm_idx);
        
        % CORRECT FOR AOM OFFSET
        sync_shots.probe_set(idx) = 2*sync_shots.wm_setpt(idx) - opts.tr.aom_freq;
        
        [~,this_tdc_idx] = closest_value(tdc_time,this_time+25);
        this_tdc_idx = this_tdc_idx ;
        sync_shots.tdc_time(idx) = tdc_time(this_tdc_idx);
        sync_shots.N_atoms(idx) = data.tdc.N_atoms(this_tdc_idx)';
    end
    
    ctrl_mask = logical(sync_shots.lv_cal)';
    % Masking for laser setpt errors
    wm_set_err = sync_shots.wm_setpt - 2*sync_shots.lv_set/1e6; % Error in MHz
    wm_set_mask = ~isnan(sync_shots.wm_setpt)';
    wm_msr_mask = (abs(wm_set_err) < opts.tr.wm_tolerance)';
    
    %Create the masks
    cal_mask = ctrl_mask;
    msr_mask = ~ctrl_mask & wm_set_mask;
    
    % Break out calibration % measurement blocks
    sync_cal = struct_mask(sync_shots,cal_mask);
    sync_msr = struct_mask(sync_shots,msr_mask);
        
    

    
    %% Write the output
    sync_data.shots = sync_shots;
    sync_data.cal = sync_cal;
    sync_data.msr = sync_msr;
    
    
    header({1,'Done.'})
    
    %% Plot diagnostics
    
    [cal_hist,cal_bin_edges]= histcounts(sync_cal.N_atoms,opts.tr.num_cal_bins);
    cal_bin_cents = 0.5*(cal_bin_edges(1:end-1)+cal_bin_edges(2:end));
    
    [wm_hist,cal_wm_bin_edges] = histcounts(wm_set_err(wm_msr_mask),opts.tr.num_cal_bins);
    wm_bin_cents = 0.5*(cal_wm_bin_edges(1:end-1)+cal_wm_bin_edges(2:end));
    
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
    
    filename1 = fullfile(opts.tr.out_dir,'timestamp_match');
    saveas(f2,[filename1,'.fig']);
    saveas(f2,[filename1,'.png'])
    
    
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
end