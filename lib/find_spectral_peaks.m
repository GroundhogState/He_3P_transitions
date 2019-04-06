function peaks = find_spectral_peaks(data,opts)
header({0,'Finding peaks...'})
    signal = 1-data.spec.signal;
    satch_mask = signal > opts.saturation_threshold;
    sig_satch = signal;
    sig_satch(satch_mask) = 1;
    freq_cut = data.spec.freq;

    df = mean(diff(data.spec.freq));
    
    smooth_out = smoothdata(sig_satch,'gaussian',opts.smooth_width);
    smooth_cut = smooth_out.*(smooth_out>opts.cutoff_thresh);    
    % Values, indices, widths, prominences
    [pks,locs,w,p] = findpeaks(smooth_cut); %
    
    
    
    %% Return output
    peaks.peak_num = (1:length(pks))';
    peaks.vals = signal(locs)';
    peaks.locs = locs';
    peaks.freqs = data.spec.freq(locs)';
    peaks.widths = (w*df)';
    peaks.prominences = p';
    
    
    %% Verbose output
    
    header({2,'%u peaks found.',length(locs)})

    
    sfigure(50123);
    clf;
    plot(data.spec.freq,signal,'g.')
    hold on
    plot(data.spec.freq,sig_satch,'b.')
    plot(freq_cut,smooth_out,'k.')
    plot(freq_cut,smooth_cut,'ro')    
    for pidx = 1:length(pks)
        loc = data.spec.freq(locs(pidx));
       plot(loc.*[1,1],[0,-0.1],'k-','LineWidth',1.0) 
       plot(loc.*[1,1]+0.5*peaks.widths(pidx)*[-1,1],[-0.05,-0.05],'k-','LineWidth',1.0) 
    end
    legend('Raw signal','forced saturation','Smoothed signal','Thresholded signal','Peak locations') 
    title('Automatic peak detection')    
header({1,'Done.'})    
    
end