function data_tdc = import_mcp_tdc(opts)
    cli_header({0,'Importing TDC data'})
    
    cache_opts=[];
    cache_opts.verbose=1;

    if ~isfield(opts.tdc,'force_load_save')
        opts.tdc.force_load_save=false;
    end
    cache_opts.force_cache_load=opts.tdc.force_load_save;
    opts.tdc=rmfield(opts.tdc,'force_load_save');

    if ~isfield(opts.tdc,'force_forc')
        opts.tdc.force_forc=false;
    end

    if opts.tdc.force_forc %if force_forc then have to skip the cache
        opts.tdc.force_reimport=true;
    end

    if ~isfield(opts.tdc,'force_reimport')
        opts.tdc.force_reimport=false;
    end
    cache_opts.force_recalc=opts.tdc.force_reimport;
    opts.tdc=rmfield(opts.tdc,'force_reimport');

    data_tdc=simple_function_cache(cache_opts,@mcp_tdc_import_core,{opts});
%     data_tdc=outputs{1};
    
    
    if isfield(opts.tdc,'mcp_post_fun')
        data_tdc = opts.tdc.mcp_post_fun(data_tdc,opts);
    else
        data_tdc = mcp_post_fun(data_tdc,opts);
    end
    if opts.tdc.plots
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
        
        
        filename = fullfile(opts.tdc.out_dir,sprintf('%s_log',mfilename));
        saveas(f,[filename,'.fig']);
        saveas(f,[filename,'.png']);
        
        
    end
    cli_header({1,'Done.'})
end


function data_tdc = mcp_post_fun(data_tdc,opts)

    % Don't need anything else yet...
    data_tdc.N_atoms = data_tdc.num_counts';

end


function mcp_tdc_data = mcp_tdc_import_core(opts)
    opts.tdc.shot_num=find_data_files(opts.tdc);
    %opts.tdc.shot_num= opts.tdc.shot_num(1:10); %debuging
    out.dir = opts.tdc.dir;
    %start up the diary of stdout
    diary([out.dir,'anal.txt'])
    %import the data
    [mcp_tdc_data,opts.tdc]=import_mcp_tdc_data(opts.tdc);

end

