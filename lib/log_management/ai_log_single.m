function ai_log_single_out=ai_log_single(args_single)
% single mode interogated over the entire range
%in
%  args_single.dir
%  args_single.fname
%  args_single.sfp.num_checks
%  args_single.sfp.thresh_cmp_peak
%  args_single.sfp.peak_dist_min_pass
%  args_single.pd.time_start
%  args_single.pd.time_stop

% outs
% ai_log_single_out.pd.mean
% ai_log_single_out.pd.std
% ai_log_single_out.pd.median
% ai_log_out.ok.sfp



%%load the data
path=strcat(args_single.dir,args_single.fname);
fid = fopen(path,'r');
raw_line = fread(fid,Inf,'*char')'; % this is the fastest string read method 3.5s for 20 files
fclose(fid);
ai_log_single_out=jsondecode(raw_line);
    
end