function elogs = check_error_logs(data,opts)
    % Accepts data after time sync and outputs a set of error check masks, including a master
    % mask. By default, masks nothing, but flags errors.
    header({0,'Checking the logs for errors...'}) 
    
    if ~isfield(opts,'pd_offset'), opts.pd_offset = -0.06; end
    if ~isfield(opts,'pd_delay'), opts.pd_delay = 0; end
    if ~isfield(opts,'pd_setpoint'), opts.pd_setpoint = .07; end
    if ~isfield(opts,'pd_setpoint_tol'), opts.pd_setpoint_tol = .03; end
    if ~isfield(opts,'pd_noise_tol'), opts.pd_noise_tol = 1e-4; end
    if ~isfield(opts,'min_atom_num'), opts.min_atom_num = 5e3; end
    
    if isfield(data,'sync')
        if isfield(data.sync,'shots')
            s_data = data.sync.shots;
            num_shots = length(s_data.tdc_time);
        else 
            warning('No time-synced data found, wtf?')
        end
    else
        warning('No time-synced data found, wtf?')
    end
    
    
    num_msr = length(data.sync.msr.tdc_time);
	wm_set = zeros(num_shots,1);
    wm_act = zeros(num_shots,1);
    wm_blue = ones(num_shots,1);
    pd_set = zeros(num_shots,1);
    pd_clean = zeros(num_shots,1);

    
    

    
    if isfield(data,'lv')
       if isfield(data.lv,'setpoint')
           wm_set = data.lv.setpoint/1e6; %Converting to MHz
       else
           warning('No setpoint recorded in LV logs!')
       end
    end
    
    if isfield(data,'wm')
       if isfield(data.wm.feedback,'actual') 
        	wm_act = data.wm.feedback.actual;
       else
           warning('Laser operating frequencies not reported!!!')
       end
       if isfield(data.wm.feedback,'blue_freq')
%            wm_blue = data.wm.feedback.blue_freq;
       end
    else
        warning('No wavemeter logs found!!!')
    end
    
    if isfield(data,'ai')     
       % Retrieve photodiode setpoint
%         pd_set = mean(window_of_ai_record);
%         pd = var(window_of_ai_record);
       % Compute noise on photodiode
       if size(data.ai.pd_data,1) > 3e4, data.ai.pd_data = data.ai.pd_data(1:3e4,:); end
       num_ai = size(data.ai.pd_data,2);
       for nai = 1:num_ai
            pd_record = (data.ai.pd_data(:,nai) - opts.pd_offset);
            pd_high = pd_record > .02;
            pd_on_val = pd_record(pd_high)';
            pd_mean(nai) = mean(pd_on_val);
            pd_var(nai) = var(pd_on_val);
       end
       pd_set = abs(pd_mean - opts.pd_setpoint) < opts.pd_setpoint_tol;
       pd_clean = pd_var < opts.pd_noise_tol;
    else
        warning('NO analog input logs found, output may not be trustworthy!')
        elogs = data.sync;
    end
    
    opts.wm_err_tolerance = 1; %MHz
    
    
       % Check atom number in calibration shots
       
       % If a measurement shot is between two dead calibration shots, ignore it
       
       
       % Check the wavemeter set point error - note calibrations have NaN setpoints!
%     elogs.wm_set_check = wm_set - wm_act < opts.wm_err_tolerance;
       % Check that the doubler is locked
    elogs.wm_blue_check = wm_blue > 0; %might need to be more sophisticated here
       % Check that the setpoint is reached
    elogs.pd_set = pd_set & pd_clean;
       % Check the intensity noise
    elogs.master = elogs.pd_set' & elogs.wm_blue_check;
    sprintf('Discarding %u/%u measurement shots in error log checks',num_msr-sum(elogs.master),num_msr)
    header({1,'Done.'})
end