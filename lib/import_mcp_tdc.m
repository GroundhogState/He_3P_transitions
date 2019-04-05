function data_tdc = import_mcp_tdc(opts_tdc)
    fwtext('Importing TDC data')
    
    data_tdc = mcp_tdc_import_core(opts_tdc);
    
%     if isfield(opts_tdc,mcp_post_fun)
        data_tdc = mcp_post_fun(data_tdc,opts_tdc);
%     end
    if opts_tdc.plots
        f = sfigure(400);
        
        subplot(2,2,1)
        plot(data_tdc.shot_num,data_tdc.num_counts,'.')
        xlabel('Shot number')
        ylabel('Number of atoms')
        title('Hit count trend')
        
        subplot(2,2,2)
        plot(data_tdc.shot_num,data_tdc.time_create_write(:,2),'.')
        
        subplot(2,1,2)
        plot(data_tdc.time_create_write(:,2),data_tdc.num_counts,'.')
        
        suptitle('DLD import diagnostics')
        
        
        filename = fullfile(opts_tdc.out_dir,sprintf('%s_log',mfilename));
        saveas(f,[filename,'.fig']);
        saveas(f,[filename,'.png']);
        
        
    end
end


function data_tdc = mcp_post_fun(data_tdc,opts_tdc)

    % Don't need anything else yet...
    data_tdc.N_atoms = data_tdc.num_counts';

end


function mcp_tdc_data = mcp_tdc_import_core(opts_tdc)

    header({0,'Importing...'})
    opts_tdc.shot_num=find_data_files(opts_tdc);
    %opts_tdc.shot_num= opts_tdc.shot_num(1:10); %debuging
    out.dir = opts_tdc.dir;
    %start up the diary of stdout
    diary([out.dir,'anal.txt'])
    %import the data
    [mcp_tdc_data,import_opts]=import_mcp_tdc_data(opts_tdc);
    header({0,'Done.'})
end

