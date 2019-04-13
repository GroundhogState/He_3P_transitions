function opts = dir_config(opts,data_dir)
    opts.dir = data_dir;
    %set up an output dir %https://gist.github.com/ferryzhou/2269380
    if (exist([opts.dir,'out'], 'dir') == 0), mkdir([opts.dir,'out']); end
    %make a subfolder with the ISO timestamp for that date
    opts.out_dir=sprintf('%sout\\%s\\',...
        opts.dir,datestr(datetime('now'),'yyyymmddTHHMMSS'));
    if (exist(opts.out_dir, 'dir') == 0), mkdir(opts.out_dir); end
    opts.lv.dir = opts.dir;
    opts.lv.out_dir = opts.out_dir;
    opts.ai.dir=opts.dir;
    opts.ai.out_dir = opts.out_dir;
    opts.wm.dir=opts.dir;
    opts.wm.out_dir = opts.out_dir;
    opts.tdc.dir = opts.dir;
    opts.tdc.out_dir = opts.out_dir;
    opts.ai.cache_single.mock_working_dir=opts.dir;
    opts.wm.cache_import.mock_working_dir=opts.dir;
    wm_logs=dir([opts.wm.dir,opts.wm.wm_log_name,'*.txt']);
    opts.wm.names={wm_logs.name};
    

    %Invoke local config
end