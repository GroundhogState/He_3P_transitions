function sync_data = match_timestamps(data,opts)

    %% Matching timestamps
    % Should:   Truncate all logs by time to fit within smallest logged
    % window
    % Plot all diagnostics
    % Extract as function/consider generalizing
    % What's with the 10MHz offset?!
    
    header({0,'Correlating timestamps...'})
    
    tdc_time = data.tdc.time_create_write(:,2);
    if ~isfield(data.ai,'error')
        ai_time = data.ai.timestamp;
    end
    lv_time = data.lv.time;
    wm_time = data.wm.feedback.posix_time;
    file_idxs = 1:length(tdc_time);
    
    
    start_time = min(lv_time);
    sync_data.start_time = start_time;
    all_setpts = unique(data.lv.setpoint);
    lv_set_range = [min(all_setpts)/1e6,max(all_setpts)/1e6]; %MHz
    mid_setpt = mean(lv_set_range); 
    
    if isfield(opts,'ignorefiles')
        if ~isnan(opts.ignorefiles)
           warning('Manually ignoring %u of %u files!!! Check opts.filemask',length(opts.ignorefiles),length(file_idxs))
           file_idxs = setdiff(file_idxs,opts.ignorefiles);
           num_files = length(file_idxs); %filemask is a list of shots to drop
        else
            num_files = length(file_idxs);
        end
    else
       num_files = length(file_idxs);
    end
    
    blue_rec = false;
    if isfield(data.wm,'blue_freq')
        blue_rec = true;
    else
       warning('Blue wavelength not recorded by WS8!!')
    end
    

%     if isfield(data.ai,'error')
% %         data.ai.ai_mask = ones(size(tdc_time));
%         warning('AI check failure - bypassing mask, results may not be trustworthy!')
%     end
        sync_shots = [];
    for idx = 1:num_files
        this_tdc_idx = file_idxs(idx); % To allow for manual masking
       % Fun job: Write fn/ modify closest_value so it can be vectorized
        this_time = tdc_time(this_tdc_idx);
        
        sync_shots.tdc_time(idx) = tdc_time(this_tdc_idx);
        sync_shots.N_atoms(idx) = data.tdc.N_atoms(this_tdc_idx)';
        
        
        [~,this_lv_idx] = closest_value(lv_time,this_time-25);
        sync_shots.lv_cal(idx) = data.lv.calibration(this_lv_idx);
        sync_shots.lv_set(idx) = data.lv.setpoint(this_lv_idx);
        sync_shots.class(idx) = data.lv.shot_class(this_lv_idx);        
        sync_shots.lv_time(idx) = this_time;

        if ~isfield(data.ai,'error')
            ai_offset_time = this_time;
            ai_diffs = ai_time-ai_offset_time;
            ai_past = ai_diffs<0;
            this_ai_idx = find(ai_past,1,'last');
            this_ai_time = ai_time(this_ai_idx);
            sync_shots.mean_power(idx) = data.ai.high_sum(this_ai_idx)/data.ai.high_time(this_ai_idx);
            sync_shots.exposure(idx) = data.ai.high_time(this_ai_idx);
            this_probe_window = data.ai.probe_window(:,this_ai_idx);
            
            [this_wm_time,this_wm_idx] = closest_value(wm_time,this_time-21);
            sync_shots.wm_time(idx) = this_wm_time;
            [up_time,up_idx] = closest_value(wm_time,this_probe_window(1));
            [down_time,down_idx] = closest_value(wm_time,this_probe_window(2));
            wm_active_log = data.wm.feedback.actual(up_idx:down_idx);
            sync_shots.wm_std(idx) = std(wm_active_log);
            sync_shots.wm_setpt(idx) = mean(wm_active_log);
        else
            sync_shots.mean_power(idx) = 1;
            sync_shots.exposure(idx) = 1;
%             this_probe_window = data.ai.probe_window(:,this_ai_idx);
            
            [this_wm_time,this_wm_idx] = closest_value(wm_time,this_time-21);
            sync_shots.wm_time(idx) = this_wm_time;
            [~,up_idx] = closest_value(wm_time,this_wm_time);
            [~,down_idx] = closest_value(wm_time,this_wm_time+.25);
            wm_active_log = data.wm.feedback.actual(up_idx:down_idx);
            sync_shots.wm_std(idx) = std(wm_active_log);
            sync_shots.wm_setpt(idx) = mean(wm_active_log);
        end
%         probe_window(idx,:) = this_probe_window;
        
        

%         if blue_rec
%             sync_shots.wm_blue(idx) = data.wm.blue_freq.value(this_wm_idx);
%         end
        % CORRECT FOR AOM OFFSET
        
%         [this_ai_time,this_ai_idx] = closest_value(ai_time,this_time-21);

%         sync_shots.ai_check(idx) = data.ai.ai_mask(this_ai_idx);
        
        sync_shots.probe_set(idx) = 2*sync_shots.wm_setpt(idx) - opts.aom_freq;
        
    end

    
    
    %%     Splitting calibration and measurements 
%     ai_mask = sync_shots.ai_check';
    ctrl_mask = logical(sync_shots.lv_cal)';

    %% Masking for laser setpt errors
%     sync_shots.wm_set_err = sync_shots.wm_setpt - sync_shots.lv_set/1e6; % Error in MHz
    sync_shots.wm_set_err = sync_shots.lv_set/1e6 - sync_shots.wm_setpt;
    wm_set_mask = abs(sync_shots.wm_set_err) < opts.check.wm_tolerance;
    wm_msr_mask = ~isnan(sync_shots.wm_set_err);
%     wm_msr_mask = (abs(sync_shots.wm_set_err(wm_set_mask)-wm_err_mean) < opts.wm_tolerance)';
    if sum(wm_msr_mask) >0
        fprintf('Wavemeter set exceeds %.2f MHz in %u measurement shots\n',opts.check.wm_tolerance, length(ctrl_mask)-sum(wm_msr_mask))
    end
  
    % Performing safety checks
    elogs = [];
%     elogs = check_error_logs(data,opts);
    

    
    %Create the masks
    cal_mask = ctrl_mask;
    if isfield(elogs,'master')
        msr_mask = ~ctrl_mask & wm_set_mask' & elogs.master;
    else
        msr_mask = ~ctrl_mask & wm_set_mask' ;
    end
        
    % Break out calibration % measurement blocks
    sync_cal = struct_mask(sync_shots,cal_mask);
    sync_msr = struct_mask(sync_shots,msr_mask);
        

    %% Write the output
    sync_data.shots = sync_shots;
    sync_data.cal = sync_cal;
    sync_data.msr = sync_msr;
    


    
    header({1,'Done.'})
    
    %% Plot diagnostics
    
    [cal_hist,cal_bin_edges]= histcounts(sync_cal.N_atoms,opts.check.num_cal_bins);
    cal_bin_cents = 0.5*(cal_bin_edges(1:end-1)+cal_bin_edges(2:end));
    
    [wm_hist,cal_wm_bin_edges] = histcounts(sync_shots.wm_set_err(wm_msr_mask),opts.check.num_cal_bins);
    wm_bin_cents = 0.5*(cal_wm_bin_edges(1:end-1)+cal_wm_bin_edges(2:end));
    
    f2=sfigure(500);
    clf
%     subplot(2,3,1)
%     area(cal_bin_cents,cal_hist)
%     xlabel('Value [V]')
%     title(sprintf('Atom number in calibration shots, N=%u',sum(cal_mask)))

    subplot(2,3,1)
    plot(sync_shots.wm_setpt,'x')
    xlabel('Shot number')
    ylabel('Measured WM freq')
    title('Probe setpt')
    
%     subplot(2,3,3)
%     area(wm_bin_cents,wm_hist)
%     xlabel('Value [MHz]')
%     title(sprintf('Laser setpoint error, N=%u',length(sync_shots.wm_set_err(wm_set_mask))))
    
%     subplot(2,3,5)
%     plot(sync_shots.lv_time-start_time,sync_shots.lv_set/1e6-data.wm.feedback.actual,'x');
%     hold on
%     plot(wm_time-start_time,data.wm.feedback.actual,'x');
%     legend('2*LV set point','WM blue setpt')  
%     xlabel('Time elapsed')
%     title('Raw channel inputs')
    
    
    subplot(2,3,2)
    plot(wm_time-start_time,'.')
    hold on
    plot(tdc_time-start_time,'.')
    plot(lv_time-start_time,'.')
    if ~isfield(data.ai,'error')
        plot(ai_time-start_time,'.')
    end
    xlabel('Shot number')
    ylabel('Time')
    title('Recorded timestamps')
    legend('WM','TDC','LV','AI')
    
    filename1 = fullfile(opts.out_dir,'timestamp_match');
    saveas(f2,[filename1,'.fig']);
    saveas(f2,[filename1,'.png'])
    
    subplot(2,3,3)
    plot(sync_shots.lv_time-start_time,sync_shots.lv_set/1e6,'x');
    hold on
    plot(wm_time-start_time,data.wm.feedback.actual,'x');
    legend('2*LV set point','WM blue setpt')  
    xlabel('Time elapsed')
    title('Raw channel inputs')
    subplot(2,3,4)
    plot(sync_shots.lv_time-start_time,sync_shots.lv_set/1e6,'x');
    hold on
    plot(sync_shots.lv_time-start_time,sync_shots.wm_setpt,'x');
    legend('2*LV set point','WM blue setpt')  
    xlabel('Time elapsed')
    ylabel(sprintf('Frequency - %3f',mid_setpt))
    legend('LabView log set','WM record')
    title('Time-matched setpoints')
    
    subplot(2,3,5)
    plot(sync_msr.lv_time,sync_msr.lv_set/1e6 - sync_msr.wm_setpt,'x')
    xlabel('Time elapsed')
    ylabel('2*LV set - WM set')
    title('Setpoint error, trimmed')
    
    subplot(2,3,6)
    plot(sync_msr.lv_time,sync_msr.probe_set,'x')
    
    title('AOM-corrected scanned range, trimmed')
    
    suptitle('Probe beam setpoint matching')
end