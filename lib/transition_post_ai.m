function post_ai = transition_post_ai(opts_ai,data_ai)

%     opts_ai.pd_offset = -.065
    
    dt = 1./opts_ai.srate;

    post_ai = data_ai;
    ai_data = cellfun(@(x) x.Data, data_ai.data,'UniformOutput',false);
    pd_data= cell2mat(cellfun(@(x) x(3,:), ai_data,'UniformOutput',false)')';
    opts_ai.pd_offset = mean(pd_data(1:1000,:));
    pd_reset= pd_data- opts_ai.pd_offset;
    pd_smooth= smoothdata(pd_reset); %Smooth out noise
    high_mask = pd_smooth > opts_ai.high_thresh;
    
    
    post_ai.high_sum = sum(high_mask.*pd_smooth);
    post_ai.high_time = dt*sum(high_mask);
    post_ai.high_mean = post_ai.high_sum./sum(high_mask);
    
    
    post_ai.pd_range = range(pd_smooth);
    iso_to_posix = @(x) posixtime(datetime(x,'InputFormat','yyyyMMdd''T''HHmmss'));
    post_ai.timestamp = cell2mat(nucellf(@(x) iso_to_posix(x.ISO_time_start_aq),data_ai.data));
%     post_ai.t0 = min(post_ai.timestamp);

    up_idx = zeros(size(high_mask,2),2);
    for shot = 1:length(up_idx)
       start_idx = find(high_mask(:,shot),1);
        end_idx = find(high_mask(:,shot),1,'last');
       if ~isempty(start_idx)
           up_idx(shot,:) = [start_idx,end_idx];
       end 
        
    end
    up_window = dt*up_idx';
    
    sample_duration = dt*size(data_ai.data{1}.Data,2);
    post_ai.probe_window = post_ai.timestamp - sample_duration + up_window;
    
    
    
%     post_ai.ai_mask = ones(size(post_ai.high_time));
%     on_mask = post_ai.high_time > opts_ai.exposure;
%     power_mask = high_mean > opts_ai.pd_setpoint;
%     post_ai.ai_mask = on_mask & power_mask;
%     fprintf('PD performance failure in %u of %u shots\n',sum(post_ai.ai_mask),length(post_ai.ai_mask));
    
    ts_offset = post_ai.timestamp(1);

    

    f=stfig('Analog import diagnostics');
    clf;

    
    
    subplot(2,2,1)
    plot(post_ai.timestamp-ts_offset,'.')
    xlabel('File number')
    ylabel('Time elapsed from start')
    suptitle('Analog import diagnostics')

    subplot(2,2,2)
    plot(post_ai.pd_range,'.')
    title('Smoothed PD range')
    xlabel('Shot number')
    ylabel('Voltage range')
    
    subplot(2,2,3)
    plot(post_ai.high_time,'.')
    title('PD on time (s)')
    
    subplot(2,2,4)
    plot(post_ai.high_mean,'.')
    title('Avg high power')

    filename = fullfile(opts_ai.out_dir,sprintf('%s_log',mfilename));
    saveas(f,[filename,'.fig']);
    saveas(f,[filename,'.png']);
        
    
end