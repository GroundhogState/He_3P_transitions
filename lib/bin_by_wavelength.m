function data = bin_by_wavelength(data,opts)
num_cats = numel(data.cat);
for cidx = 1:num_cats
cli_header({0,'Binning data for presentation'})
    sync_msr = data.cat{cidx}.data;
    all_setpts = unique(sync_msr.probe_set);
    all_setpts = all_setpts(~isnan(all_setpts));
    min_set = min(sync_msr.probe_set);
    max_set = max(sync_msr.probe_set);
    
    
    num_freq_bins = round(range(all_setpts/1e6));
    freq_bin_size = 1;%mean(diff(sort_setpt));
    fbin_edges = linspace(min_set-freq_bin_size,max_set+freq_bin_size,num_freq_bins+1);%/1e6 %MHz
    if isfield(opts,'num_freq_bins')
        if ~isnan(opts.num_freq_bins) 
            num_freq_bins = opts.num_freq_bins;
        end
    end
    if isfield(opts,'freq_bin_size') %overrides number-of-bins
        if ~isnan(opts.freq_bin_size) 
            freq_bin_size = opts.freq_bin_size;
            fbin_edges = min_set-freq_bin_size:freq_bin_size:max_set+freq_bin_size;%/1e6 %MHz
            num_freq_bins = length(fbin_edges)-1;
        end
    end
    
    
    [freq_sorted,freq_idx] = sort(sync_msr.probe_set);
    shots_sort = struct_mask(sync_msr,freq_idx);
    %test this
    
    
    freq_stats.freq=nan(numel(fbin_edges)-1,1);
%      1     9    39    52    56
    for ii=1:numel(fbin_edges)-1
       %fast_sorted_mask(freq_sorted,fbin_edges(ii),fbin_edges(ii+1))
       fmask_temp = shots_sort.probe_set > fbin_edges(ii) & shots_sort.probe_set < fbin_edges(ii+1);
       if sum(fmask_temp) >0
%            iis = msr_set_id(fmask_temp);
           shots_temp = struct_mask(shots_sort,fmask_temp);
           %test this
           
           freq_stats.freq(ii) = nanmean(shots_temp.probe_set);
           
           freq_stats.num_shots(ii) = sum(fmask_temp);
           freq_stats.freq_std(ii) = nanstd(shots_temp.probe_set);
           freq_stats.sig_cal(ii) = nanmean(shots_temp.signal);
           freq_stats.sig_std(ii) = nanstd(shots_temp.signal);
           
       end
    end
%     [msr_setpts ,msr_set_id]= sort(sync_msr.probe_set);  
    
%     for ii=1:num_freq_bins
%        fmask_temp = msr_setpts > fbin_edges(ii) & msr_setpts < fbin_edges(ii+1);
%        if sum(fmask_temp) >0
%            iis = msr_set_id(fmask_temp);
%            shots_temp = struct_mask(sync_msr,iis);
%            freq_stats.num_shots(ii) = sum(fmask_temp);
%            freq_stats.freq(ii) = nanmean(shots_temp.probe_set);
%            freq_stats.freq_std(ii) = nanstd(shots_temp.probe_set);
%            freq_stats.sig_cal(ii) = nanmean(shots_temp.signal);
%            freq_stats.sig_std(ii) = nanstd(shots_temp.signal);
%            
%        end
%     end
    
    
    
    if isfield(freq_stats,'sig_err')
        freq_stats = rmfield(freq_stats,'sig_err');
        freq_stats = rmfield(freq_stats,'freq_err');
    end
    freq_stat_mask = freq_stats.freq > 0;
    freq_stats = struct_mask(freq_stats,freq_stat_mask);
    freq_stats.freq_err = freq_stats.freq_std./sqrt(freq_stats.num_shots);
    freq_stats.sig_err = freq_stats.sig_std./sqrt(freq_stats.num_shots);
    
   
    data.cat{cidx}.freq_stats = freq_stats;
    cli_header({1,'Done.'})
    
    X = freq_stats.freq;
    Y = freq_stats.sig_cal;
    X_err =  freq_stats.freq_err;
    Y_err = freq_stats.sig_err;
    
    f1=sfigure(2347+cidx);
    subplot(2,1,1)
%     plot(
%     plot(data.cat{cidx}.freq_stats.freq,data.cat{cidx}.freq_stats.sig_cal,'ko-')
    plot(sync_msr.probe_set,sync_msr.signal,'k.-')
    xlabel('Frequency (MHz)')
    ylabel('N/N_0')
    title('Raw spectral data')
%     errorbar(
    subplot(2,1,2)
    errorbar(X,Y,Y_err,Y_err,X_err,X_err)
    title(sprintf('Spectrum with %.1fMHz bins',freq_bin_size))
    xlabel('Frequency (MHz)')
    ylabel('N/N_0')
    
    filename1 = fullfile(opts.out_dir,'binned_spectrum');
    saveas(f1,[filename1,'.fig']);
    saveas(f1,[filename1,'.png'])
end
end