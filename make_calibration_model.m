function calib  = make_calibration_model(data,opts)    
    
    
    calib = interp1(data.tr.sync.cal.tdc_time,data.tr.sync.cal.N_atoms,data.tr.sync.msr.tdc_time);      

    % Plotting the calibration model
        
    plot_pred = opts.pred_freq.*1e3-data.tr.fits.lorz_fit(1);
    plot_raw_T = data.tr.sync.msr.tdc_time-data.tr.start_time;
    plot_raw_X = data.tr.fits.freqs_centred;
    plot_raw_Y = data.tr.sync.msr.N_atoms;
    
    
    plot_cal_T = data.tr.sync.cal.tdc_time-data.tr.start_time;
    plot_cal_X = data.tr.sync.cal.probe_set;
    plot_cal_Y = data.tr.sync.cal.N_atoms;
    
    plot_mdl_Y = data.tr.sync.msr.N_atoms./calib';
    
    
    plot_fit_X = linspace(min([data.tr.fits.freqs_centred,plot_pred-1]),max([data.tr.fits.freqs_centred,plot_pred+1]),1e4);
    plot_fit_Y_lor = data.tr.fits.lor_mdlfun(data.tr.fits.lor_cen_prms,plot_fit_X);
    plot_fit_Y_gaus = data.tr.fits.gaus_mdlfun(data.tr.fits.gaus_cen_prms,plot_fit_X);
    
    f1=sfigure(5000);
    clf;
    subplot(4,1,1)
    plot(plot_raw_T,plot_raw_Y,'.')
    hold on
    plot(plot_cal_T,plot_cal_Y,'.')
    plot(plot_raw_T,calib)
    xlabel('Time elapsed')
    ylabel('Counts')
    legend('Measurement shots','calibration shots','Model')
    title('Raw data')

    subplot(4,1,2)
    plot(plot_raw_T,plot_mdl_Y,'.')
    xlabel('Time elapsed')
    ylabel('N ratio')
    title('Model-calibrated signal')    
    suptitle('Calibration model')
    
    subplot(4,1,3)
    plot(data.tr.sync.shots.lv_time-data.tr.start_time,data.tr.sync.shots.wm_setpt,'.')
    xlabel('Time elapsed')
    ylabel('Raw WM set point')
    
    
    subplot(4,1,4)
    plot(plot_raw_X,plot_mdl_Y,'.')
    hold on
    plot(plot_fit_X,plot_fit_Y_lor,'r-')
    plot(plot_fit_X,plot_fit_Y_gaus,'k--')
    xlabel(sprintf('f-%3.5f [MHz]'))
    ylabel('N ratio')
    ylim([-0.1,1.1])
    title('Model-calibrated signal')    
    suptitle('Calibration model')
    legend('Calibrated data','Lorentzian fit','Gaussian fit')
    
    
    filename2 = fullfile(opts.out_dir,sprintf('%s_diagnostic',mfilename));
    saveas(f1,[filename2,'.fig']);
    saveas(f1,[filename2,'.png'])


    end