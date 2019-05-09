function sync_data = match_timestamps(data,opts)

    %% Matching timestamps
    % Should:   Truncate all logs by time to fit within smallest logged
    % window
    % Plot all diagnostics
    % Extract as function/consider generalizing
    % What's with the 10MHz offset?!
    
    cli_header({0,'Correlating timestamps...'})
    
    tdc_time = data.tdc.time_create_write(:,1);
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
    
    if ~isfield(opts.ai,'force_idx_match'), opts.ai.force_idx_match=false; end
    if opts.ai.force_idx_match
        warning('Overriding AI timestamp matching!!')
    end
    if ~isfield(data.ai,'error')
       ai_time = data.ai.timestamp;
       ai_stash = data.ai; 
    else
       ai_stash.timestamp = tdc_time;
    end


        sync_shots = [];
    for idx = 1:num_files
        if ~isempty(ai_stash.timestamp)
            this_tdc_idx = file_idxs(idx); % To allow for manual masking
            % Fun job: Write fn/ modify closest_value so it can be vectorized
            this_time = tdc_time(this_tdc_idx);
            
            [this_ai_time,this_ai_idx] = closest_value(ai_stash.timestamp,this_time-21);
            if opts.ai.force_idx_match         
                this_ai_idx = idx;
            end
            if abs(this_ai_time - (this_time-21)) > 10 % files too far apart
                disp('Analog input too far from nearest shot')
                continue
            else

                sync_shots.tdc_time(idx) = tdc_time(this_tdc_idx);
                sync_shots.N_atoms(idx) = data.tdc.N_atoms(this_tdc_idx)';     

                [this_lv_time,this_lv_idx] = closest_value(lv_time,this_time-20);
                sync_shots.lv_cal(idx) = data.lv.calibration(this_lv_idx);
                sync_shots.lv_set(idx) = data.lv.setpoint(this_lv_idx);
                sync_shots.class(idx) = data.lv.shot_class(this_lv_idx);        
                sync_shots.lv_time(idx) = this_lv_time;
                [this_wm_time,this_wm_idx] = closest_value(wm_time,this_lv_time+2);
                sync_shots.wm_time(idx) = this_wm_time;
                    if ~isfield(data.ai,'error')
% 

%                             this_ai_idx = find(ai_past,1,'last');
%                         end
                            sync_shots.mean_power(idx) = ai_stash.high_mean(this_ai_idx);
                            sync_shots.exposure(idx) = ai_stash.high_time(this_ai_idx);
                            this_probe_window = ai_stash.probe_window(:,this_ai_idx);

                            [up_time,up_idx] = closest_value(wm_time,this_lv_time+3); 
%                             if length(this_probe_window) < 2 || diff(this_probe_window) == 0
                                down_idx = up_idx + 1;
%                             else
%                                 [down_time,down_idx] = closest_value(wm_time,this_probe_window(2));
%                             end
                            wm_active_log = data.wm.feedback.actual(up_idx:down_idx);
                            sync_shots.wm_std(idx) = std(wm_active_log);
                            sync_shots.wm_setpt(idx) = mean(wm_active_log);
                            sync_shots.idxs(idx,:) = [this_tdc_idx,this_lv_idx,this_wm_idx,this_ai_idx];
                            sync_shots.times(idx,:) = [this_time,this_lv_time,this_wm_time,this_ai_time]-start_time;
                            sync_shots.ai_time(idx) = this_ai_time;
                            % Remove this analog file from record so it can't be reused
                            if ~opts.ai.force_idx_match
                                ai_filt = zeros(1,length(ai_stash.timestamp));
                                ai_filt(this_ai_idx) = 1;
                                ai_stash = struct_mask(ai_stash,~logical(ai_filt));
                            end
                    else
                        sync_shots.mean_power(idx) = 1;
                        sync_shots.exposure(idx) = 1;
            %             this_probe_window = data.ai.probe_window(:,this_ai_idx);
                        sync_shots.idxs(idx,:) = [this_tdc_idx,this_lv_idx,this_wm_idx];
                        [~,up_idx] = closest_value(wm_time,this_wm_time);
                        [~,down_idx] = closest_value(wm_time,this_wm_time+.25);
                        sync_shots.times(idx,:) = [this_time,this_lv_time,this_wm_time]-start_time;
                    end
%                 [this_wm_time,this_wm_idx] = closest_value(wm_time,this_time-21);
%                 sync_shots.wm_time(idx) = this_wm_time;
                wm_active_log = data.wm.feedback.actual(up_idx:down_idx);
                sync_shots.wm_trace{idx} = wm_active_log;
                sync_shots.wm_std(idx) = std(wm_active_log);
                sync_shots.wm_setpt(idx) = mean(wm_active_log);
                sync_shots.probe_set(idx) = 2*sync_shots.wm_setpt(idx) - opts.aom_freq;
            end 
        end
    end

    
    
    %%     Splitting calibration and measurements 
%     ai_mask = sync_shots.ai_check';
    ctrl_mask = logical(sync_shots.lv_cal)';

    %% Masking for laser setpt errors
%     sync_shots.wm_set_err = sync_shots.wm_setpt - sync_shots.lv_set/1e6; % Error in MHz
    sync_shots.wm_set_err = sync_shots.lv_set/1e6 - sync_shots.wm_setpt;
    wm_set_mask = abs(sync_shots.wm_set_err) < opts.check.wm_tolerance;
    wm_msr_mask = ~isnan(sync_shots.wm_set_err);
 

    
    %Create the masks
    cal_mask = ctrl_mask;
    msr_mask = ~ctrl_mask & wm_set_mask' ;

        
    % Break out calibration % measurement blocks
    sync_cal = struct_mask(sync_shots,cal_mask);
    sync_msr = struct_mask(sync_shots,msr_mask);

    wm_err_ok = abs(sync_msr.wm_set_err) < opts.check.wm_tolerance;
    fprintf('Wavemeter set exceeds %.2f MHz in %u/%u measurement shots\n',opts.check.wm_tolerance, sum(~wm_err_ok) ,length(wm_err_ok))

  

    %% Write the output
    sync_data.shots = sync_shots;
    sync_data.cal = sync_cal;
    sync_data.msr = sync_msr;
    


    
    cli_header({1,'Done.'})
    
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
%     plot(wm_time-start_time,'.')
    plot(sync_shots.tdc_time-start_time,'o')
    hold on
    plot(sync_shots.lv_time-start_time,'.')
    if ~isfield(data.ai,'error')
        plot(sync_shots.ai_time-start_time,'*')
    end
    xlabel('Shot number')
    ylabel('Time')
    title('Recorded timestamps')
    legend('TDC','LV','AI')
    
    subplot(2,3,3)
    plot(sync_shots.lv_time-start_time,sync_shots.lv_set/1e6,'ro');
    hold on
    plot(wm_time-start_time,data.wm.feedback.actual,'.');
    legend('2*LV set point','WM blue setpt')  
    xlabel('Time elapsed')
    title('Raw channel inputs')
    
    subplot(2,3,4)
    plot(sync_shots.lv_time-start_time,sync_shots.wm_setpt,'.');
    hold on
    plot(sync_shots.lv_time-start_time,sync_shots.lv_set/1e6,'ko');
    legend('WM blue setpt','2*LV set point')  
    xlabel('Time elapsed')
    ylabel(sprintf('Frequency - %3f',mid_setpt))
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
    
        filename1 = fullfile(opts.out_dir,'timestamp_match');
    saveas(f2,[filename1,'.fig']);
    saveas(f2,[filename1,'.png'])
end