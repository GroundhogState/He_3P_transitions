function post_ai = transition_post_ai(opts_ai,data_ai)

    post_ai = data_ai;
    ai_data = cellfun(@(x) x.Data, data_ai.data,'UniformOutput',false);
    post_ai.pd_data= cell2mat(cellfun(@(x) x(3,:), ai_data,'UniformOutput',false)')';
    post_ai.pd_range = range(post_ai.pd_data);
    iso_to_posix = @(x) posixtime(datetime(x,'InputFormat','yyyyMMdd''T''HHmmss'));
    post_ai.timestamp = cell2mat(nucellf(@(x) iso_to_posix(x.ISO_time_start_aq),data_ai.data));
    post_ai.t0 = min(post_ai.timestamp);

end