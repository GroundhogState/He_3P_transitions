function data = bin_by_wavelength(data,opts)
num_cats = numel(data.cat);
for cidx = 1:num_cats
header({0,'Binning data for presentation'})
    sync_msr = data.cat{cidx}.data;
    
    all_setpts = unique(data.lv.setpoint);
    min_set = min(sync_msr.probe_set);%2*min(all_setpts);
    max_set = max(sync_msr.probe_set);%2*max(all_setpts);
    
    num_freq_bins = round(range(all_setpts/1e6));
    if isfield(opts,'num_freq_bins')
        if ~isnan(opts.num_freq_bins) 
            num_freq_bins = opts.num_freq_bins;
        end
    end
    freq_gap = 0;%mean(diff(sort_setpt));
    fbin_edges = linspace(min_set-freq_gap,max_set+freq_gap,num_freq_bins+1);%/1e6 %MHz
    [msr_setpts ,msr_set_id]= sort(sync_msr.probe_set);

    for ii=1:num_freq_bins
       fmask_temp = msr_setpts > fbin_edges(ii) & msr_setpts < fbin_edges(ii+1);
       if sum(fmask_temp) >0
           iis = msr_set_id(fmask_temp);
           shots_temp = struct_mask(sync_msr,iis);
           freq_stats.num_shots(ii) = sum(fmask_temp);
           freq_stats.freq(ii) = nanmean(shots_temp.probe_set);
           freq_stats.freq_std(ii) = nanstd(shots_temp.probe_set);
           freq_stats.sig(ii) = nanmean(shots_temp.N_atoms);
           freq_stats.sig_cal(ii) = nanmean(shots_temp.N_atoms./shots_temp.calib);
           freq_stats.sig_std(ii) = nanstd(shots_temp.N_atoms./shots_temp.calib);
           
       end
    end
    freq_stat_mask = abs(freq_stats.freq) > 0;
    freq_stats = struct_mask(freq_stats,freq_stat_mask);
    freq_stats.freq_err = freq_stats.freq_std./sqrt(freq_stats.num_shots);
    freq_stats.sig_err = freq_stats.sig_std./sqrt(freq_stats.num_shots);
    
   
    data.cat{cidx}.freq_stats = struct_mask(freq_stats,logical(freq_stat_mask));
    header({1,'Done.'})
    
    X = freq_stats.freq;
    Y = freq_stats.sig_cal;
    X_err =  freq_stats.freq_err;
    Y_err = freq_stats.sig_err;
    
    sfigure(2347+cidx);
%     plot(data.cat{cidx}.freq_stats.freq,data.cat{cidx}.freq_stats.sig_cal,'ko-')
%     errorbar(
    errorbar(X,Y,Y_err,Y_err,X_err,X_err)
    suptitle('Spectrum with 1MHz bins')
    xlabel('Frequency (MHz)')
    ylabel('N/N_0')
end
end