
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
addpath(genpath(core_folder));
zeeman_folder = fullfile(fileparts(this_folder),'Zeeman_splitting\');
addpath(genpath(zeeman_folder));

hebec_constants
% initialize variables
% Note that master config also calls a local override if it exists



root_dir = "\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\2019_He_Spectroscopy\spectroscopy_data\final_scans\";
data_dirs = {"20190417_5^1D2_mj_1_itc_both_qwp_146_overnight"
            "20190425_5^3D1_both_ITC_1MHz_step"
            "20190429_5^3D_2_3_qwp_142"
            "20190429_5^3D_2_3_qwp_142_tight_scan_centre_peak"
            "20190429_5^3S_1_qwp_146_both_stage"};
data_paths =  cellfun(@(x) fullfile(root_dir,x),data_dirs);
out_file = "output_and_options.mat";
n_dirs = numel(data_dirs);
for p=2:5
    dpath = data_paths{p};
    opts = master_transition_config(dpath);
    opts.fidx = p;
    opts.font = 'Times';
    if exist(lopt_file,'file')==2
        addpath(dpath)
        opts = local_opts(opts,dpath);
        rmpath(dpath)
    end
    old_data = get_last_output(dpath,out_file);
    lopt_file = fullfile(dpath,'local_opts.m');
    data.cat = old_data.out_data.data;
    opts.out_dir = "C:\Users\jacob\Documents\Projects\He_3P_transitions\figs\";
    combo_plots(data,opts);
end



%%
function loaded_data=get_last_output(varargin)
   dpath = varargin{1};
    out_file=varargin{2};
    outpath = fullfile(dpath,'out');
    files = dir(outpath);
    fnames = {files.name};
    is_subdir = ~ismember(fnames,{'.','..'});
    fnames = fnames(is_subdir);
    fstamps = cell2mat(cellfun(@(x) str2num(strrep(x,'T','')),fnames,'uni',0)');
    [~,f_order] = sort(fstamps,'descend');
    n_subdirs = numel(fnames);
    %Starting from latest output, see whether output file was saved
    for pp = 1:n_subdirs
        dir_idx = f_order(pp);
        this_name = fnames{dir_idx};
        this_dir = fullfile(outpath,this_name);
        target_file = fullfile(this_dir,out_file);
        if exist(target_file,'file')
            fprintf('File FOUND in %s\n',this_dir)
            loaded_data = load(target_file);
            fprintf('Data loaded from %s\n',this_dir)
            break
        else
        end 
    end
end