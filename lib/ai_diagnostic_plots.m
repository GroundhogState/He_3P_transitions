function ai_diagnostic_plots(out_ai,opts_ai)

    ts_offset = out_ai.timestamp(1);

    f=sfigure(100);
    clf

    subplot(1,2,1)
    plot(out_ai.timestamp-ts_offset,'.')
    xlabel('File number')
    ylabel('Time elapsed from start')
    suptitle('Analog import diagnostics')

    subplot(1,2,2)
    plot(out_ai.pd_range,'.')
    title('PD data')
    xlabel('Sample idx')
    ylabel('Voltage range')

    filename = fullfile(opts_ai.out_dir,sprintf('%s_log',mfilename));
    saveas(f,[filename,'.fig']);
    saveas(f,[filename,'.png']);
        
end