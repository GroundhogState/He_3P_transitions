function data_tdc = import_mcp_tdc(opts_tdc)
    header({0,'Importing TDC data'})
    
    cache_opts=[];
    cache_opts.verbose=1;

    if ~isfield(opts_tdc,'force_load_save')
        opts_tdc.force_load_save=false;
    end
    cache_opts.force_cache_load=opts_tdc.force_load_save;
    opts_tdc=rmfield(opts_tdc,'force_load_save');

    if ~isfield(opts_tdc,'force_forc')
        opts_tdc.force_forc=false;
    end

    % if ~isfield(opts_tdc,'no_save')
    %     %to be completed, requires modification to function_cache
    % end

    if opts_tdc.force_forc %if force_forc then have to skip the cache
        opts_tdc.force_reimport=true;
    end

    if ~isfield(opts_tdc,'force_reimport')
        opts_tdc.force_reimport=false;
    end
    cache_opts.force_recalc=opts_tdc.force_reimport;
    opts_tdc=rmfield(opts_tdc,'force_reimport');

    outputs=function_cache(cache_opts,@mcp_tdc_import_core,{opts_tdc});
    data_tdc=outputs{1};
    
%     data_tdc = mcp_tdc_import_core(opts_tdc);
    
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
    header({1,'Done.'})
end


function data_tdc = mcp_post_fun(data_tdc,opts_tdc)

    % Don't need anything else yet...
    data_tdc.N_atoms = data_tdc.num_counts';

end


function mcp_tdc_data = mcp_tdc_import_core(opts_tdc)
    opts_tdc.shot_num=find_data_files(opts_tdc);
    %opts_tdc.shot_num= opts_tdc.shot_num(1:10); %debuging
    out.dir = opts_tdc.dir;
    %start up the diary of stdout
    diary([out.dir,'anal.txt'])
    %import the data
    [mcp_tdc_data,opts_tdc]=import_mcp_tdc_data(opts_tdc);

end

