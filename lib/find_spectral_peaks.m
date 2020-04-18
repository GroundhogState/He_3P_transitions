function peaks = find_spectral_peaks(spec,opts)
cli_header(0,'Finding peaks...');
    signal = spec.signal;
%     satch_mask = signal < opts.peak.saturation_threshold;
%     sig_satch = signal;
%     sig_satch(satch_mask) = 1;
    freq_cut = spec.freq;
    sig_satch = signal;

    df = mean(diff(spec.freq));
    
    smooth_out = smoothdata(sig_satch,'gaussian',opts.peak.smooth_width);
    smooth_cut = smooth_out.*(smooth_out>opts.peak.cutoff_thresh);    
    % Values, indices, widths, prominences
    [pks,locs,w,p] = findpeaks(smooth_cut); %
    
    
    
    %% Return output
    peaks.peak_num = (1:length(pks))';
    peaks.vals = signal(locs)';
    peaks.locs = locs';
    peaks.freqs = spec.freq(locs)';
    peaks.widths = (w*df)';
    peaks.prominences = p';
    
    
    %% Verbose output
    
    cli_header(2,'%u peaks found.',length(locs));
    fnum = 50123;
    if isfield(opts,'fig_idx')
        fnum = fnum + 10*opts.fig_idx;
    end
    f1=stfig('Peak detection');
    subplot(2,1,opts.fig_idx);
    cla;
    plot(spec.freq,signal,'k.:')
    hold on
%     plot(spec.freq,sig_satch,'g.')
    plot(freq_cut,smooth_out,'b*-')
    plot(freq_cut,smooth_cut,'rx')    
    for pidx = 1:length(pks)
        loc = spec.freq(locs(pidx));
        val = peaks.vals(pidx);
        plot(loc.*[1,1],[val,-0],'k-','LineWidth',1.0) 
        plot(loc.*[1,1]+0.5*peaks.widths(pidx)*[-1,1],0.5*val*[1,1],'k-','LineWidth',1.0) 
    end
%     xtickformat('%u')
    xlabel('Frequency (MHz)')
    legend('Raw signal','Smoothed signal','Thresholded signal','Peak locations') 
    title(sprintf('Detecting peaks in stage %u',opts.fig_idx))
    suptitle('Automatic peak detection')    
    
    filename1 = fullfile(opts.out_dir,'peak_detection');
    saveas(f1,[filename1,'.fig']);
    saveas(f1,[filename1,'.png'])
cli_header(1,'Done.');    
    
end