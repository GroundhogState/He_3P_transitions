function out_ai=ai_log_import(opts_ai)

fwtext('IMPORTING DATA')

if opts_ai.verbose > 0
    header({0,'Importing analog logs...'})
end
%a simple wrapper for the below ai_log_import that uses the matlab function cache
cache_opts=opts_ai.cache_import;

%limit the scope but retain the structure
data_sub = [];
out_ai=simple_function_cache(cache_opts,@ai_log_import_core,{opts_ai,data_sub});

% Kind of violates the concept of this being a 'pure' import script, so could compute this
% elsewhere. It's easy to start computing things in here as diagnostic which end up getting used
% later... But for the sake of component testing, this is here. 
% Could be abstracted further  “¯\_(?)_/¯“

if isfield(opts_ai,'post_fun')
    header({2,'Post-processing...'})
    out_ai = opts_ai.post_fun(opts_ai,out_ai);
    header({2,'Done'})
    
    % Can be slow to plot large imports, but a) not likely to be enabled often and b) only called once
    % per import
    if opts_ai.plots % Add clause so this automatically disabled if reloading from cache?        
        ai_diagnostic_plots(out_ai)
    end
end

if opts_ai.verbose > 0
    header({1,'Done!'})
end

end

function ai_log_out=ai_log_import_core(opts_ai,data)
%ai_log_import - imports analog input log and checks if the probe beam pd signal is ok and that the laser is single
%mode by measuring the distance (in pzt voltage) between the scanning FP pd peaks
%results are placed into a convenient strucure with mat file cache
%will load data from a cashed version if opts_ai has not changed
%the output is a well aranged structure to be added into the data structure
%the data structure is not included in the save to prevent double saving of the data
%at the end of the import or load cashe the approapriate feilds are added to data
%Only checks if sm if pd is ok/all previous checks of sm are ok
%
% Syntax:  [data,import_opts]=import_data(import_opts)
%
    
% Known BUGS/ Possible Improvements
%   -need more dam speed!
%       -reading and jsondecode are the main problems the rest is very fast
%   - vector for setpt
%   -Improved documentation
%   - Replace file counter with progress bar to reduce CLI write calls
% Author: Bryce Henson
% email: Bryce.Henson[a circle]live.com  %YOU MUST INCLUDE '[matlab][ai_log_import]' in the subject line OR I WILL NOT REPLY
% Last revision:2018-09-30

%estimate of the sfp scan time,used to set the window and the smoothing

args_single=opts_ai.args_single;

%------------- BEGIN CODE --------------

dir_read=dir([opts_ai.dir,opts_ai.log_name,'*.txt']);
ai_log_out.file_names={dir_read.name};

%set up for the ai_log_single 
args_single.dir=opts_ai.dir;

cache_opts = opts_ai.cache_single;


iimax=size(ai_log_out.file_names,2); %the number of ai logs that have been identified
% Manual setting for number of files to import 
if ~isnan(opts_ai.num_files)
    iimax = min(iimax,opts_ai.num_files);
end
%initalize outputs
%loop over all the ai_logs

if opts_ai.cache_single_import
    fprintf('Caching all imports, may be slow.\n')
else
    fprintf('Caching only after import.\n')
end

fprintf('processing ai log for files %04u:%04u',iimax,0)

for ii=1:iimax
%     try
    if mod(ii,10) == 0
        fprintf('\b\b\b\b%04i',ii)
    end
    if ii==iimax || ii==1
        cache_opts.clean_cache=true; %clean cache at start and end
    else
        cache_opts.clean_cache=false;
    end
    cache_opts.clean_cache=false;
    fname=ai_log_out.file_names{ii};
    args_single.fname=fname;
    
%   Caching the import can take a long time if there are many files to import. The time to save the
%   cache outputs is massive, so you may want to disable this option and just cache the entire
%   input instead of single calls. 

    if opts_ai.cache_single_import
        ai_log_out.data{ii}=simple_function_cache(cache_opts,@ai_log_single,{args_single});    
    else
        ai_log_out.data{ii}=ai_log_single(args_single);
    end

end %loop over files
fprintf('\nDone\n')

end%function


