function peaks = find_spectral_peaks(data,opts)

   [freq_sorted,freq_order] = sort(data.sync.msr.probe_set);
    sig_sorted = 1-data.sync.msr.N_atoms./data.sync.msr.calib';
    sig_sorted = sig_sorted(freq_order);
    satch_mask = sig_sorted > opts.saturation_threshold;
    sig_satch = sig_sorted;
    sig_satch(satch_mask) = 1;
    freq_cut = freq_sorted;
    
    smooth_out = smoothdata(sig_satch,'gaussian',opts.smooth_width);
    smooth_cut = smooth_out.*(smooth_out>opts.cutoff_thresh);    
    % Values, indices, widths, prominences
    [pks,locs,w,p] = findpeaks(smooth_cut); %
    
    
    %% Return output
    peaks.vals = pks;
    peaks.locs = locs;
    peaks.freqs = freq_sorted(locs);
    peaks.widths = w;
    peaks.prominences = p;
    
    
    %% Verbose output
    
    fprintf('%u peaks found:\n',length(locs))
    for pidx = 1:length(pks)
        fprintf('Peak %u: %.3f MHz\n',pidx,peaks.freqs(pidx))
    end
    
    sfigure(50123);
    clf;
    plot(freq_sorted,sig_sorted,'g.')
    hold on
    plot(freq_sorted,sig_satch,'b.')
    plot(freq_cut,smooth_out,'k.-')
    plot(freq_cut,smooth_cut,'r.')    
    for pidx = 1:length(pks)
       plot(freq_sorted(locs(pidx)).*[1,1],[0,-0.1],'b-','LineWidth',1.0) 
    end
    legend('Raw signal','forced saturation','Smoothed signal','Thresholded signal','Peak locations') 
    suptitle('Automatic peak detection')    
    
    
end