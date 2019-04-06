function freq_stats = bin_by_wavelength(data,opts)
header({0,'Binning data for presentation'})
    sync_msr = data.tr.sync.msr;
    
    
    all_setpts = unique(data.lv.setpoint);
    min_set = min(sync_msr.probe_set);%2*min(all_setpts);
    max_set = max(sync_msr.probe_set);%2*max(all_setpts);
    
    num_freq_bins = round(range(all_setpts/1e6));
    if ~isfield(opts,'num_freq_bins')
        if ~isnan(opts.num_freq_bins) 
            num_freq_bins = opts.num_freq_bins;
        end
    end
    freq_gap = 0;%mean(diff(sort_setpt));
    fbin_edges = linspace(min_set-freq_gap,max_set+freq_gap,num_freq_bins+1);%/1e6 %MHz
    [msr_setpts ,msr_set_id]= sort(sync_msr.probe_set);
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
           freq_stats.freq(ii) = nanmean(shots_temp.probe_set);
           freq_stats.sig(ii) = nanmean(shots_temp.N_atoms);
           freq_stats.freq_err(ii) = nanstd(shots_temp.probe_set);
           freq_stats.num_shots(ii) = sum(fmask_temp);
           freq_stats.sig_cal(ii) = nanmean(shots_temp.N_atoms./shots_temp.calib');
           freq_stats.sig_std(ii) = nanstd(shots_temp.N_atoms./shots_temp.calib');
           freq_stats.sig_err(ii) = freq_stats.sig_std(ii)/sqrt(freq_stats.num_shots(ii));
           freq_stat_mask(ii) = abs(freq_stats.freq(ii)) > 0;
       end
    end
    freq_stats = struct_mask(freq_stats,logical(freq_stat_mask));
    header({1,'Done.'})
end