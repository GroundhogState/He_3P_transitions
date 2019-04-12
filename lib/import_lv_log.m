function data = import_lv_log(opts)
    header({0,'Importing LV log...'})
    lv_log=[];
    lv_log.dir = strcat(opts.lv.dir,'log_LabviewMatlab.txt');
    fid = fopen(lv_log.dir);
    lv_log.cell=textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    lv_log.cell=lv_log.cell{1};
    for ii=1:size(lv_log.cell,1)
        if ~isequal(lv_log.cell{ii},'') %catch the empty case
            if contains(lv_log.cell{ii},'measure_probe')
                line_cells=textscan(lv_log.cell{ii},'%f %s %s %s %f %s %u','Delimiter',',');
                lv_log.setpoints(ii)=line_cells{5};
                lv_log.probe_calibration(ii)=false;
                lv_log.iter_nums(ii)=line_cells{7};
                lv_log.shot_class{ii} = erase(line_cells{4}{1},'measure_probe_');
            elseif contains(lv_log.cell{ii},'calibrate')
                line_cells=textscan(lv_log.cell{ii},'%f %s %s %s %s %u','Delimiter',',');
                lv_log.setpoints(ii)=NaN;
                lv_log.probe_calibration(ii)=true;
                lv_log.iter_nums(ii)=line_cells{6};
                lv_log.shot_class{ii} = 'calibration';
            else %deals with the legacy case (only 20180813_CW_AL_tuneout_scan)
                line_cells=textscan(lv_log.cell{ii},'%f %s %s %s %f %s %u','Delimiter',',');
                lv_log.setpoints(ii)=line_cells{5};
                lv_log.probe_calibration(ii)=false;
                lv_log.iter_nums(ii)=line_cells{7};
                lv_log.shot_class{ii} = 'null';
            end
            lv_log.posix_times(ii)=line_cells{1};
            lv_log.iso_times{ii}=line_cells{2};
        end
    end
    data.setpoint=lv_log.setpoints*1e6; %convert to hz
    data.time=lv_log.posix_times;
    data.shot_num=lv_log.iter_nums;
    data.calibration=lv_log.probe_calibration;
    data.shot_class = lv_log.shot_class;
    header({1,'Done!'})
    
    if opts.lv.plots
        t0 = min(data.time);
       
        f =  sfigure(200);
        clf;
        subplot(2,2,1)
        plot(data.shot_num,data.time-t0,'.')
        xlabel('Shot number')
        ylabel('Elapsed time')
        ylabel('Shot time')
        subplot(2,2,2)
        plot(data.shot_num,data.setpoint,'.')
        xlabel('Shot number')
        ylabel('WM setpoint')
        title('Wavelength scan')
        subplot(2,2,[3 4])
        plot(data.shot_num,data.calibration,'.')
        xlabel('Shot number')
        ylabel('Calibration?')
        title('Calibration mask')

        suptitle('LabView log diagnostics')
        
        filename = fullfile(opts.lv.out_dir,sprintf('%s_log',mfilename));
        saveas(f,[filename,'.fig']);
        saveas(f,[filename,'.png']);
    end
end