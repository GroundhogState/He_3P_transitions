function wm=wm_log_import(opts)

cli_header(0,'Importing WM data');

cache_opts = opts.wm.cache_import;
cache_opts.verbose = 1;
opts.wm.force_recalc = false;
wm=simple_function_cache(cache_opts,@wm_log_import_core,{opts});


if opts.wm.plots
    wm_time = wm.feedback.posix_time;
    t0 = min(wm_time);
    wm_set = wm.feedback.setpt;
    wm_act = wm.feedback.actual;
    
    f = stfig('Wavemeter import diagnostics');
    clf
    subplot(2,1,1)
    plot(wm_time-t0,wm_set,'k.')
    hold on
    plot(wm_time-t0,wm_act,'r.')
    legend('Set','Actual')
    xlabel('Time elapsed (s)')
    ylabel('WM setpoint (MHz)')
    title('Blue freq setpoint')
    
    subplot(2,1,2)
    plot(wm_time-t0,wm_set-wm_act,'k.')
    title('WM error')
    xlabel('Time elapsed')
    ylabel('Set-actual freq (MHz)')
    
    suptitle('Wavemeter diagnostics')
    
    
    filename = fullfile(opts.wm.out_dir,sprintf('%s_log',mfilename));
    saveas(f,[filename,'.fig']);
    saveas(f,[filename,'.png']);
end

cli_header(1,'Done');

end

function wm_log=wm_log_import_core(opts)
%read the json formated wm-laser logfile
%the json is formated as {"posix_time":num,"iso_time":"2018-08-25T20:07:01.763","oper":{operation a structure}}
%want to ingest into a structure
%wm_log.operationa.posix_time
%wm_log.operationa.

%injest these logs into a 
%struct with each type of event
%wm_log.feedback
%with wm_log.feedback.time wm_log.feedback.set_wav
wm_log=struct();
add_to_struct=true;
fprintf('Importing %u wavemeter-laser feedback log files',size(opts.wm.names,2))
iimax=size(opts.wm.names,2);
ndigits = floor(log10(iimax));
% zero_string = '';
% for p=1:ndigits
%     zero_string = [zero_string,'0']
% end
% fprintf(['Importing file ',zero_string,'/',iimax,'\n'])
for ii=1:iimax
    path=strcat(opts.wm.dir,opts.wm.names{ii});
    fid = fopen(path,'r');
    wm_log_file_cells=textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
%     for p=1:ndigits+1
%         fprintf('\b\b')
%     end
    fprintf(['%u/%u\n',ii,iimax])
    fprintf('\nFile %03u/%03u importing %03uk lines:%03uk',ii,iimax,round(size(wm_log_file_cells{1},1)*1e-3),0)
    %now process these lines and place the entries into a feild of the wm_log struct depending on the operation performed
    jjmax=size(wm_log_file_cells{1},1);
    if ~isnan(opts.wm.num_logs)
        jjmax = min(jjmax, opts.wm.num_logs);
    end
    for jj=1:jjmax
%         if mod(jj,1e4)==0,fprintf('\b\b\b\b%03uk',round(jj*1e-3)),end  %fprintf('\b\b\b\b%04u',jj)
        if ~isempty(wm_log_file_cells{1}{jj})
            try
                line_tmp=wm_log_file_cells{1}{jj};
                
                if jj==jjmax %deal with the last line of the log which specifies the path of the new log, with the incorrect json syntax (single slash instead of double)
                    line_tmp=strrep(line_tmp,'\\','\');%change the new correct format(double slash) to the old incorrect single slash
                    line_tmp=strrep(line_tmp,'\','\\');%change the new correct format(double slash) to the old incorrect single slash
                end
                formated_line=jsondecode(line_tmp);
            catch
                warning('json encode faiiled file %u line %u \n\n',ii,jj)
                fprintf('\nFile %03u/%03u importing %03uk lines:%03uk',ii,iimax,round(size(wm_log_file_cells{1},1)*1e-3),0)
                add_to_struct=false;
            end

            %iso_to_posix=posixtime(datetime(formated_line.iso_time,'InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSS'));
            %if iso_to_posix~=formated_line.posix_time
                %warning('times do not match taking ISO one')
                %formated_line.posix_time=iso_to_posix;
            %end
            if add_to_struct
                %find the type of operation
                layer1_feilds=fieldnames(formated_line);
                fmask=~(strcmp(fieldnames(formated_line),'posix_time') | strcmp(fieldnames(formated_line),'iso_time'));
                if sum(fmask)>1
                    error('more than one operation')
                end
                operation=layer1_feilds(fmask);
                operation=operation{1};
                if isstruct(formated_line.(operation)) %deal with my bad formating of having operations have values directly
                    layer2_feilds=fieldnames(formated_line.(operation)); %loop over all the feild names in the struct result
                    if sum(strcmp(operation,fieldnames(wm_log)))==0 %initalize this operation feild
                        %fprintf('make feild\n')
                        wm_log.(operation)=[];
                        for fn=layer2_feilds'
                            wm_log.(operation).(fn{1})={};
                        end
                        wm_log.(operation).posix_time=[];
                    end
                    %append to each matrix in the structure
                    mat_idx=size(wm_log.(operation).posix_time,2)+1;
                    wm_log.(operation).posix_time(mat_idx)=formated_line.posix_time;
                    for fn=layer2_feilds'
                        if isstruct(formated_line.(operation).(fn{1}))
                            for fn_layer3=fieldnames(formated_line.(operation).(fn{1}))'
                                wm_log.(operation).(fn_layer3{1}){mat_idx}=formated_line.(operation).(fn{1}).(fn_layer3{1});
                            end
                        else
                            wm_log.(operation).(fn{1}){mat_idx}=formated_line.(operation).(fn{1});
                        end
                    end
                else
                    if sum(strcmp(operation,fieldnames(wm_log)))==0 %initlialize
                        wm_log.(operation)=[];
                        wm_log.(operation).posix_time=[];
                        wm_log.(operation).value=[];
                    end
                    mat_idx=size(wm_log.(operation).posix_time,2)+1;
                    wm_log.(operation).value{mat_idx}=formated_line.(operation);
                    wm_log.(operation).posix_time(mat_idx)=formated_line.posix_time;
                end
            else
                add_to_struct=true;
            end
        end
    end
end
fprintf('\nCleaning up output structure...\n')
wm_log=clean_log_structure(wm_log,[]); %itteratively defined cleaner



end


%log example
% {"posix_time":1535227621.763,"iso_time":"2018-08-25T20:07:01.763","feedback":{"setpt":362858686.395,"actual":362858685.698,"Res":49.817107,"Int":49.817803,"slew_lim":0}}
% {"posix_time":1535227621.834,"iso_time":"2018-08-25T20:07:01.834","feedback":{"setpt":362858686.395,"actual":362858685.766,"Res":49.817942,"Int":49.818571,"slew_lim":0}}
% {"posix_time":1535227621.895,"iso_time":"2018-08-25T20:07:01.895","feedback":{"setpt":362858686.395,"actual":362858685.989,"Res":49.818741,"Int":49.819147,"slew_lim":0}}
% {"posix_time":1535227621.970,"iso_time":"2018-08-25T20:07:01.970","feedback":{"setpt":362858686.395,"actual":362858686.006,"Res":49.819232,"Int":49.819621,"slew_lim":0}}
% {"posix_time":1535227622.042,"iso_time":"2018-08-25T20:07:02.042","feedback":{"setpt":362858686.395,"actual":362858685.741,"Res":49.819948,"Int":49.820602,"slew_lim":0}}
% {"posix_time":1535227622.101,"iso_time":"2018-08-25T20:07:02.101","feedback":{"setpt":362858686.395,"actual":362858686.185,"Res":49.820690,"Int":49.820900,"slew_lim":0}}
% {"posix_time":1535227622.159,"iso_time":"2018-08-25T20:07:02.159","feedback":{"setpt":362858686.395,"actual":362858686.135,"Res":49.820952,"Int":49.821211,"slew_lim":0}}
% {"posix_time":1535227622.211,"iso_time":"2018-08-25T20:07:02.211","feedback":{"setpt":362858686.395,"actual":362858686.263,"Res":49.821232,"Int":49.821364,"slew_lim":0}}
% {"posix_time":1535227622.285,"iso_time":"2018-08-25T20:07:02.285","feedback":{"setpt":362858686.395,"actual":362858686.198,"Res":49.821372,"Int":49.821569,"slew_lim":0}}
% {"posix_time":1535227622.344,"iso_time":"2018-08-25T20:07:02.344","feedback":{"setpt":362858686.395,"actual":362858686.406,"Res":49.821564,"Int":49.821553,"slew_lim":0}}
% {"posix_time":1.535227622431E+9,"iso_time":"2018-08-25T20:07:02.431","get_status":{"status":0,"wavelength":826.537842,"temperature":221.641769,"temperature_status":"off","etalon_lock":"on","etalon_voltage":42.221703,"cavity_lock":"off","resonator_voltage":98.018082,"ecd_lock":"on","ecd_voltage":98.018082,"output_monitor":1.015346,"etalon_pd_dc":0.832644,"dither":"on"}}
% {"posix_time":1.535227622431E+9,"iso_time":"2018-08-25T20:07:02.431","read_all_adc":{"status":0,"channelCount":18,"channel0":"Ref.Cav.thermistor","value0":1.991442,"units0":"V","channel1":"PCB thermistor","value1":0.782691,"units1":"V","channel2":"Aux in","value2":0.00155,"units2":"V","channel3":"TP8","value3":0.0001,"units3":"V","channel4":"Etalon PD DC","value4":0.832644,"units4":"V","channel5":"Ref.Cav.PD","value5":0.000111,"units5":"V","channel6":"Output PD","value6":1.015346,"units6":"V","channel7":"Input PD","value7":0.000149,"units7":"V","channel8":"ECD PD2","value8":0.530835,"units8":"V","channel9":"ECD PD1","value9":0.592374,"units9":"V","channel10":"ECD output","value10":0.839705,"units10":"V","channel11":"ECD tuning","value11":92.036438,"units11":"V","channel12":"Ref.Cav.ext.tuning","value12":98.861794,"units12":"V","channel13":"Etalon tuning","value13":42.221703,"units13":"V","channel14":"Resonator tuning","value14":98.018082,"units14":"V","channel15":"DAC 1 ADC loopback","value15":2.208591,"units15":"V","channel16":"MONITOR A","value16":0,"units16":"V","channel17":"MONITOR B","value17":0,"units17":"V"}}
% {"posix_time":1535227622.498,"iso_time":"2018-08-25T20:07:02.498","feedback":{"setpt":362858686.395,"actual":362858686.101,"Res":49.821606,"Int":49.821900,"slew_lim":0}}
% {"posix_time":1535227622.571,"iso_time":"2018-08-25T20:07:02.571","feedback":{"setpt":362858686.395,"actual":362858685.769,"Res":49.823202,"Int":49.823828,"slew_lim":0}}
% {"posix_time":1535227622.654,"iso_time":"2018-08-25T20:07:02.654","feedback":{"setpt":362858686.395,"actual":362858686.069,"Res":49.823978,"Int":49.824303,"slew_lim":0}}
% {"posix_time":1535227622.714,"iso_time":"2018-08-25T20:07:02.714","feedback":{"setpt":362858686.395,"actual":362858686.055,"Res":49.824528,"Int":49.824868,"slew_lim":0}}
% {"posix_time":1535227622.764,"iso_time":"2018-08-25T20:07:02.764","feedback":{"setpt":362858686.395,"actual":362858686.157,"Res":49.824916,"Int":49.825154,"slew_lim":0}}
% {"posix_time":1535227622.824,"iso_time":"2018-08-25T20:07:02.824","feedback":{"setpt":362858686.395,"actual":362858686.219,"Res":49.825154,"Int":49.825330,"slew_lim":0}}
% {"posix_time":1535227622.872,"iso_time":"2018-08-25T20:07:02.872","feedback":{"setpt":362858686.395,"actual":362858686.154,"Res":49.825378,"Int":49.825618,"slew_lim":0}}
% {"posix_time":1535227622.956,"iso_time":"2018-08-25T20:07:02.956","feedback":{"setpt":362858686.395,"actual":362858686.202,"Res":49.825607,"Int":49.825799,"slew_lim":0}}
% {"posix_time":1535227623.029,"iso_time":"2018-08-25T20:07:03.029","feedback":{"setpt":362858686.395,"actual":362858686.522,"Res":49.825710,"Int":49.825583,"slew_lim":0}}
% {"posix_time":1535227623.090,"iso_time":"2018-08-25T20:07:03.090","feedback":{"setpt":362858686.395,"actual":362858686.190,"Res":49.825678,"Int":49.825882,"slew_lim":0}}
% {"posix_time":1535227623.139,"iso_time":"2018-08-25T20:07:03.139","feedback":{"setpt":362858686.395,"actual":362858686.277,"Res":49.825908,"Int":49.826026,"slew_lim":0}}
% {"posix_time":1535227623.197,"iso_time":"2018-08-25T20:07:03.197","feedback":{"setpt":362858686.395,"actual":362858686.124,"Res":49.826015,"Int":49.826286,"slew_lim":0}}
% {"posix_time":1535227623.280,"iso_time":"2018-08-25T20:07:03.280","feedback":{"setpt":362858686.395,"actual":362858685.941,"Res":49.826368,"Int":49.826822,"slew_lim":0}}
