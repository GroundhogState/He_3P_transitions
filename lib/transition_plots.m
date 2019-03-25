function plot_output = transition_plots(data,opts)    
%     pd_range_ylim = opts.pd_range_ylim;
    if ~isfield(opts,'num_cal_bins') opts.tr.num_cal_bins = 50; end
    
    % % Compute stuff

    shot_by_setpt = data.tr.shot_by_setpt;
    pd_range_cal = data.tr.pd_range_cal;
    pd_data = data.ai.pd_data;
    
    
    setpts = cellfun(@(x) x.setpoint, shot_by_setpt);
    mid_setpt = nanmedian(setpts);
    setpt_plot = setpts - mid_setpt;
    setpt_plot_MHz = setpt_plot/1e6;
    stats_plot = cell2mat(cellfun(@(x) x.stat, shot_by_setpt,'UniformOutput',false));
    pd_range_plot = stats_plot(:,1);
    pd_se_plot = stats_plot(:,2)./sqrt(stats_plot(:,3));
%     cal_bin_edges = linspace(pd_range_ylim(1),pd_range_ylim(2),opts.num_cal_bins+1);
%     cal_bin_cents = 0.5*(cal_bin_edges(1:end-1)+cal_bin_edges(2:end));

    
    
    cal_mean = stats_plot(:,1);
    cal_std = stats_plot(:,2);
    cal_se = cal_std./stats_plot(:,3);

    plt_range = [min(setpt_plot),max(setpt_plot)];
    plt_range_MHz = plt_range/1e6;


    patch_xvals = [setpt_plot_MHz',fliplr(setpt_plot_MHz)'];
    patch_yv_pos = [cal_mean',fliplr(cal_mean+cal_se)'];
    patch_yv_neg = [cal_mean',fliplr(cal_mean-cal_se)'];
    
    % % Produce plots
% if opts.tr.plot

% end
end