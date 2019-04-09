function data = auto_peak_detect(data,opts)

    num_cats = numel(data.cat);
    for cat_idx = 1:num_cats
            data_cat = data.cat{cat_idx}.data;
        % % Automatic peak detection
            [spec.freq,spec.order] = sort(data_cat.probe_set);
            sig_sorted = data_cat.N_atoms./data_cat.calib;
            spec.signal = sig_sorted(spec.order);
            opts.fig_idx = cat_idx;
            peaks = find_spectral_peaks(spec,opts);

            data.cat{cat_idx}.peaks = peaks;
            data.cat{cat_idx}.spec = spec;
            
        for pidx = 1:length(peaks.vals)
            header({1,'Peak %u data from findpeaks:',pidx})
            fprintf('Centre frequency           %.3f\n',peaks.freqs(pidx));
            fprintf('FWHM                       %.3f\n',peaks.widths(pidx));
            fprintf('Relative strength          %.2f\n',peaks.vals(pidx));
        end
    end
    
    if opts.peaks.plot
       for cat_idx = 1:num_cats
    %% Output
        midpoint = median(data.cat{cat_idx}.spec.freq);
        f=sfigure(1337+cat_idx);
        clf;
        subplot(3,1,[1 2])
        plot(data.cat{cat_idx}.spec.freq-midpoint,1-data.cat{cat_idx}.spec.signal,'k:')
        hold on
        for pidx = 1:length(data.cat{cat_idx}.peaks.vals)
           centre = data.cat{cat_idx}.peaks.freqs(pidx)-midpoint;
           height = data.cat{cat_idx}.peaks.vals(pidx);
           plot(centre*[1,1],[-0.20,height],'b-','Linewidth',2) 
           plot(centre+0.5*data.cat{cat_idx}.peaks.widths(pidx)*[-1,1],0.5*height*[1,1],'b-','Linewidth',2)
        end

        xlabel(sprintf('f - %.2f MHz',midpoint))
        ylabel('1-N_{probe}/N_{cal}')
%         suptitle([opts.tr_name,sprintf(' spectral data, ITC stage %u',cat_idx)])

        tabledata = {'Peak number','Frequency','Width','Strength'};
        for pidx = 1:length(data.cat{cat_idx}.peaks.vals)
           new_row = {pidx,sprintf('%.2f',data.cat{cat_idx}.peaks.freqs(pidx)),...
               sprintf('%.2f',data.cat{cat_idx}.peaks.widths(pidx)),sprintf('%.2f',data.cat{cat_idx}.peaks.vals(pidx))};
            tabledata = [tabledata;new_row];
        end
        uit = uitable(f,'Data',tabledata,'Position',[20 20 462 200]);

        fname=fullfile(opts.out_dir,sprintf('%s_spectrum_findpeaks',opts.tr_name));
        saveas(f,[fname,'.fig']);
        saveas(f,[fname,'.png'])


    end 
    end

end