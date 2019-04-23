function calib  = make_calibration_model(data,opts)    
    
    header({0,'Making calibration model'})
 

    % Plotting the calibration model
        
    plot_raw_T = data.sync.msr.tdc_time-data.start_time;
    plot_raw_X = data.sync.msr.probe_set;
    plot_raw_Y = data.sync.msr.N_atoms;
    
    
    plot_cal_T = data.sync.cal.tdc_time-data.start_time;
    plot_cal_X = data.sync.cal.probe_set;
    plot_cal_Y = data.sync.cal.N_atoms;
    
    calib = interp1(data.sync.cal.tdc_time,data.sync.cal.N_atoms,data.sync.msr.tdc_time);     
    plot_mdl_Y = data.sync.msr.N_atoms./calib;
    
    f1=sfigure(5000);
    clf;
    subplot(3,1,1)
    plot(plot_raw_T,plot_raw_Y,'.')
    hold on
    plot(plot_cal_T,plot_cal_Y,'.')
    plot(plot_raw_T,calib)
    xlabel('Time elapsed')
    ylabel('Counts')
    legend('Measurement shots','calibration shots','Model')
    title('Raw data')

   
    
    
    subplot(3,1,2)
    plot(plot_raw_T,plot_mdl_Y,'.')
    xlabel('Time elapsed')
    ylabel('N ratio')
    title('Calibrated time series')    
%     suptitle('Calibration model')
    
    
    
    
    subplot(3,1,3)
    sfigure(5001)
    plot(plot_raw_X,plot_mdl_Y,'.')
    hold on
    
    xlabel(sprintf('f (MHz)'))
    ylabel('N ratio')
    ylim([-0.1,1.1])
    title('Calibrated spectra')    
%     suptitle('Calibration model')
%     
%     filename2 = fullfile(opts.out_dir,sprintf('%s_diagnostic',mfilename));
%     saveas(f1,[filename2,'.fig']);
%     saveas(f1,[filename2,'.png'])

    header({1,'Done.'})

    end