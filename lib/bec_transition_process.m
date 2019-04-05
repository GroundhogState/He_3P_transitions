function data = bec_transition_process(data,opts_tr)


    %% %%%%%%%%%%%% PRECONDITIONING THE DATA 
%     num_files = numel(data.tdc.shot_num);
    tdc_time = data.tdc.time_create_write(:,2);
%     ai_time = data.ai.timestamp;
    lv_time = data.lv.time;
    wm_time = data.wm.feedback.posix_time;
    
    all_setpts = unique(data.lv.setpoint);
    lv_set_range = [min(all_setpts)/1e6,max(all_setpts)/1e6]; %MHz
    mid_setpt = mean(lv_set_range); %MHz
    
    data.tr.start_time = min(lv_time);
    num_files = length(lv_time);

    %% %%%%%%%%% COMPUTATION ON THE DATA
    
    % Match the timestamps    
    data.tr.sync = match_timestamps(data,opts_tr);

    % Create a simple calibration model
    data.tr.sync.msr.calib = make_calibration_model(data,opts_tr);
   
    % % Grouping by wavelength (Independent variable)
    data.tr.freq_stats = bin_by_wavelength(data,opts_tr);
  
    
    %% Fitting the data
   
    % Finding the peaks
    opts_tr.cutoff_thresh = 0.1;
    opts_tr.smooth_width = 20;
    opts_tr.saturation_threshold = 0.975;
    data.tr.peaks = find_spectral_peaks(data.tr,opts_tr);

    % Fit the peaks
    data.tr.fits = fit_peaks(data,opts_tr);


if opts_tr.plot
    %% Plotting
    % % Compute things
    plot_pred = opts_tr.pred_freq.*1e3-data.tr.fits.gaus_fit(1);
    
    plot_set_X = 2*sort(all_setpts);
    plot_set_X = plot_set_X(~isnan(plot_set_X))/1e6; %MHz

    plot_fit_X = linspace(min([data.tr.fits.freqs_centred,plot_pred-1]),max([data.tr.fits.freqs_centred,plot_pred+1]),1e4);
    plot_fit_Y_lor = data.tr.fits.lor_mdlfun(data.tr.fits.lor_cen_prms,plot_fit_X);
    plot_fit_Y_gaus = data.tr.fits.gaus_mdlfun(data.tr.fits.gaus_cen_prms,plot_fit_X);
        
    plot_sig_X = data.tr.freq_stats.freq-data.tr.fits.lorz_fit(1);
    plot_sig_Y = data.tr.freq_stats.sig_cal;
    plot_sig_Y_err = data.tr.freq_stats.sig_err;
    plot_sig_Y_err = plot_sig_Y_err.*~isnan(plot_sig_Y_err);
    
    % Plotting the result
    f3=sfigure(501);
    clf;
    errorbar(plot_sig_X,plot_sig_Y,plot_sig_Y_err,'.')
    hold on
    plot(plot_fit_X,plot_fit_Y_lor,'r-')
    plot(plot_fit_X,plot_fit_Y_gaus,'k--')
    plot(plot_pred.*[1,1],[0,2],'b-','LineWidth',2.0)
    xlim(1.1*[min([plot_sig_X;plot_pred-1]),max([plot_sig_X;plot_pred+1])])
    ylim([0,1.15])
    xlabel(sprintf('f - %6f (MHz)',data.tr.fits.lorz_fit(1)))
    ylabel('Atom Number Ratio')
    legend('Frequency-binned data','Lorentzian fit','Gaussian fit','Theory Value','Location','Best')

    suptitle([opts_tr.tr_name,sprintf(' absorption peak @ %.2f±(%.3f) MHz, FWHM~%.2fMHz',...
        data.tr.fits.lorz_fit(1),data.tr.fits.lorz_fit(2),data.tr.fits.lorz_fit(3))])
    title(sprintf('Vacuum wavelength %.10f nm ± %.6f fm, 2*%u shots',...
        data.tr.fits.lorz_fit(5),1e6*data.tr.fits.lorz_fit(6),length(data.tr.sync.msr.tdc_time)))
    
    % Saving 
    filename3 = fullfile(opts_tr.out_dir,sprintf('%s_spectrum',opts_tr.tr_name));
    saveas(f3,[filename3,'.fig']);
    saveas(f3,[filename3,'.png'])
    
end

end