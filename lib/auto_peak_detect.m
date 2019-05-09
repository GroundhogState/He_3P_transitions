function data = auto_peak_detect(data,opts)

    num_cats = numel(data.cat);
    for cat_idx = 1:num_cats
            data_cat = data.cat{cat_idx};
        % % Automatic peak detection
%             [spec.freq,spec.order] = sort(data_cat.probe_set);
%             sig_sorted = data_cat.N_atoms./data_cat.calib;
%             spec.signal = sig_sorted(spec.order);
            
%             [spec.freq,spec.order] = sort(data_cat.freq_stats.freq);
%             sig_sorted = data_cat.freq_stats.sig_cal;

            [spec.freq,spec.order] = sort(data_cat.data.probe_set);
            sig_to_sort = data_cat.data.signal;
            spec.signal = sig_to_sort(spec.order);
            
            opts.fig_idx = cat_idx;
            peaks = find_spectral_peaks(spec,opts);

            data.cat{cat_idx}.peaks = peaks;
            data.cat{cat_idx}.spec = spec;
            
        for pidx = 1:length(peaks.vals)
            cli_header({1,'Peak %u data from findpeaks:',pidx})
            fprintf('Centre frequency           %.3f\n',peaks.freqs(pidx));
            fprintf('FWHM                       %.3f\n',peaks.widths(pidx));
            fprintf('Relative strength          %.2f\n',peaks.vals(pidx));
%         end
    end
    

end