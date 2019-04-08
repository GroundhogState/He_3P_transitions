function data = present_data(data,opts)
num_cats = numel(data.cat);
for cat_idx = 1:num_cats
        data_cat = data.cat{cat_idx}.data;
    % % Automatic peak detection
        [data.spec.freq,data.spec.order] = sort(data_cat.probe_set);
        sig_sorted = data_cat.N_atoms./data_cat.calib';
        data.spec.signal = sig_sorted(data.spec.order);
        data.peaks = find_spectral_peaks(data,opts);

        data.cat{cat_idx}.peaks = data.peaks;
        data.cat{cat_idx}.spec = data.spec;
        
    %% Output
        midpoint = median(data.spec.freq);
        f=sfigure(1337+cat_idx);
        clf;
        subplot(3,1,[1 2])
        plot(data.spec.freq-midpoint,1-data.spec.signal,'k:')
        hold on
        for pidx = 1:length(data.peaks.vals)
           centre = data.peaks.freqs(pidx)-midpoint;
           height = data.peaks.vals(pidx);
           plot(centre*[1,1],[-0.20,height],'b-','Linewidth',2) 
           plot(centre+0.5*data.peaks.widths(pidx)*[-1,1],0.5*height*[1,1],'b-','Linewidth',2)
        end

        xlabel(sprintf('f - %.2f MHz',midpoint))
        ylabel('1-N_{probe}/N_{cal}')
        suptitle([opts.tr_name,sprintf(' spectral data, ITC stage %u',cat_idx)])

        tabledata = {'Peak number','Frequency','Width','Strength'};
        for pidx = 1:length(data.peaks.vals)
           new_row = {pidx,sprintf('%.2f',data.peaks.freqs(pidx)),...
               sprintf('%.2f',data.peaks.widths(pidx)),sprintf('%.2f',data.peaks.vals(pidx))};
            tabledata = [tabledata;new_row];
        end
        uit = uitable(f,'Data',tabledata,'Position',[20 20 462 200]);

        fname=fullfile(opts.out_dir,sprintf('%s_spectrum_findpeaks',opts.tr_name));
        saveas(f,[fname,'.fig']);
        saveas(f,[fname,'.png'])

        for pidx = 1:length(data.peaks.vals)
            header({1,'Peak %u data from findpeaks:',pidx})
            fprintf('Centre frequency           %.3f\n',data.peaks.freqs(pidx));
            fprintf('FWHM                       %.3f\n',data.peaks.widths(pidx));
            fprintf('Relative strength          %.2f\n',data.peaks.vals(pidx));
        end
    end
end