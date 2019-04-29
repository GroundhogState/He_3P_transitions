function calib  = make_calibration_model(data,opts)    
    
    header({0,'Making calibration model'})
    calib = data.sync.msr;

    % Plotting the calibration model
        
    plot_raw_T = data.sync.msr.tdc_time-data.sync.start_time;
    plot_raw_X = data.sync.msr.probe_set;
    plot_raw_Y = data.sync.msr.N_atoms;
    
    
    plot_cal_T = data.sync.cal.tdc_time-data.sync.start_time;
    plot_cal_X = data.sync.cal.probe_set;
    plot_cal_Y = data.sync.cal.N_atoms;
    

%     calib = interp1(data.sync.cal.tdc_time,data.sync.cal.N_atoms,data.sync.msr.tdc_time);     
%     plot_mdl_Y = data.sync.msr.N_atoms./calib;
%     
    % Cull out shots with too few counts based on interpolation of calibration
%     mdl_mask = data.sync.cal.N_atoms > opts.check.min_counts;
    calib.mdl = interp1(data.sync.cal.tdc_time,data.sync.cal.N_atoms,data.sync.msr.tdc_time);     
    calib.signal = calib.mdl-data.sync.msr.N_atoms;
    calib.low_count_mask = calib.mdl > opts.check.min_counts;
    calib.probe_set = data.sync.msr.probe_set;
%     calib.signal = calib.signal(calib.low_count_mask);

    
%     sync_data.msr.calib = data.cal.mdl;
%     sync_data.msr.diff = data.sync.msr.calib - data.sync.msr.N_atoms;
%     sync_data.msr.normdiff = 1-(data.sync.msr.calib - data.sync.msr.N_atoms)./data.sync.msr.calib;
%     sync_data.msr.ratio = data.sync.msr.N_atoms./data.sync.msr.calib;
    
    plot_mdl_Y = calib.signal(calib.low_count_mask);
    plot_mdl_X = data.sync.msr.probe_set(calib.low_count_mask);
    
    f1=sfigure(5000);
    clf;
    subplot(3,1,1)
    plot(plot_raw_T,plot_raw_Y,'.')
    hold on
    plot(plot_cal_T,plot_cal_Y,'.')
    plot(plot_raw_T,calib.mdl)
%     plot(plot_raw_T,calib.mdl)
    xlabel('Time elapsed')
    ylabel('Counts')
    legend('Measurement shots','calibration shots','Model')
    title('Raw data')

    
    subplot(3,1,2)
    plot(calib.signal,'.')
    hold on
    plot(calib.signal(~calib.low_count_mask),'rx')
    xlabel('Shot number')
    ylabel('N ratio')
    title('Calibrated time series')    
    legend('Raw data','N<min')
    
    
    subplot(3,1,3)
    plot(plot_mdl_X,plot_mdl_Y,'.')
    hold on

    xlabel(sprintf('f (MHz)'))
    ylabel('N ratio')

    title('Calibrated spectra')    

    filename2 = fullfile(opts.out_dir,sprintf('%s_diagnostic',mfilename));
    saveas(f1,[filename2,'.fig']);
    saveas(f1,[filename2,'.png'])

    header({1,'Done.'})

    end