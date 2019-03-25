function data_tdc = import_mcp_tdc(opts_tdc)
    fwtext('Importing TDC data')
    
    data_tdc = mcp_tdc_import_core(opts_tdc);
    
%     if isfield(opts_tdc,mcp_post_fun)
        data_tdc = mcp_post_fun(data_tdc,opts_tdc);
%     end
    
    
end


function data_tdc = mcp_post_fun(data_tdc,opts_tdc)

    % Don't need anything else yet...
    data_tdc.N_atoms = data_tdc.num_counts';

end


function mcp_tdc_data = mcp_tdc_import_core(opts_tdc)

    header({0,'Importing...'})
    opts_tdc.shot_num=find_data_files(opts_tdc);
    %opts_tdc.shot_num= opts_tdc.shot_num(1:10); %debuging

    %set up an output dir %https://gist.github.com/ferryzhou/2269380
    if (exist([opts_tdc.dir,'out'], 'dir') == 0), mkdir([opts_tdc.dir,'out']); end
    %make a subfolder with the ISO timestamp for that date
    out.dir=sprintf('%sout\\%s\\',...
        opts_tdc.dir,datestr(datetime('now'),'yyyymmddTHHMMSS'));
    if (exist(out.dir, 'dir') == 0), mkdir(out.dir); end
    opts.global.out_dir=out.dir;
    %start up the diary of stdout
    diary([out.dir,'anal.txt'])
    %import the data
    [mcp_tdc_data,import_opts]=import_mcp_tdc_data(opts_tdc);
    header({0,'Done.'})
end

