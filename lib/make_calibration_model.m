function calib  = make_calibration_model(data,opts)    
    
    cli_header(0,'Making calibration model');
    calib = data.sync.msr;

    % Plotting the calibration model
        
    plot_raw_T = data.sync.msr.tdc_time-data.sync.start_time;
    plot_raw_X = data.sync.msr.probe_set;
    plot_raw_Y = data.sync.msr.N_atoms;
    
    
    plot_cal_T = data.sync.cal.tdc_time-data.sync.start_time;
    plot_cal_X = data.sync.cal.probe_set;
    plot_cal_Y = data.sync.cal.N_atoms;
    

    calib.mdl = interp1(data.sync.cal.tdc_time,data.sync.cal.N_atoms,data.sync.msr.tdc_time);     
    calib.number_loss = calib.mdl-data.sync.msr.N_atoms;
    calib.signal = calib.number_loss./calib.mdl;
    calib.low_count_mask = calib.mdl > opts.check.min_counts;
    calib.probe_set = data.sync.msr.probe_set;

    plot_mdl_Y = calib.signal(calib.low_count_mask);
    plot_mdl_X = data.sync.msr.probe_set(calib.low_count_mask);
    
    f1=stfig('Calibration model');
    clf;
    subplot(2,1,1)
    plot(plot_raw_T,plot_raw_Y,'k.')
    hold on
    plot(plot_cal_T,plot_cal_Y,'b.')
    plot(plot_raw_T,calib.mdl,'r')
%     plot(plot_raw_T,calib.mdl)
    xlabel('Time elapsed')
    ylabel('Counts')
    legend('Measurement shots','calibration shots','Model')
    title('Raw data')

    
    subplot(2,1,2)
    plot(calib.signal,'k.')
    hold on
    plot(calib.signal(~calib.low_count_mask),'rx')
    xlabel('Shot number')
    ylabel('N ratio')
    title('Calibrated time series')    
%     legend('Raw data','N<min')
    
%     subplot(3,1,3)
%     plot(plot_mdl_X,plot_mdl_Y,'k.')
%     hold on
%     xlabel(sprintf('f (MHz)'))
%     ylabel('N ratio')
%     title('Calibrated spectra')    
    suptitle('Atom loss calibration')

    filename2 = fullfile(opts.out_dir,sprintf('%s_diagnostic',mfilename));
    saveas(f1,[filename2,'.fig']);
    saveas(f1,[filename2,'.png'])

    cli_header(1,'Done.');

    end