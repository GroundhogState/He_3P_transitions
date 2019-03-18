%A tutorial on how to get up and running doing analysis on our BEC data
% lets begin
% it is common practice for the highest level file to have some user options in it
% for our uses this is a very reasonable practice, for bigger projects it could get a tad out of hand
%%%%------------------------start user var-------------------------------%
%first we need to point at the data directory
import_opts.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
%then what the MCP-TDC data files are called 
% the practice after 2016 is for it to just be 'd' eg d123.txt
import_opts.file_name='d';
%all the other feilds are optional, and have some sensible defaults
%in order of commonly used
%   the spatio-temporal limits on the data import
%import_opts.txylim=[[0,2];[-30e-3, 30e-3];[-30e-3, 30e-3]];     %tight XY lims to eliminate hot spot from destroying pulse widths
%   the rotation angle in radians
%import_opts.dld_xy_rot=0.61;
%import_opts.force_reimport=false;
%import_opts.force_forc=false;
%%%%------------------------done user var-------------------------------%

%% Setting Up The Enviorment
%add all subfolders to the path
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));

%call the constants function that makes some globals
hebec_constants 

%% Importing Data
%now we can find all the data files in that directory
import_opts.shot_num=find_data_files(import_opts);
%and import them all into a structure, this function will also be cached in this directory so if you call the
%same import_opts it will just load the cache
data=import_mcp_tdc_data(import_opts);

%% Explore The Data Structure
% we have the fields
%            shot_num: [1×4140 double]
%            num_counts: [1×4140 double]
%            counts_txy: {1×4140 cell}
%            time_create_write: [4140×2 double]
% the one we care the most about is counts_txy
% we can look at the first file's distribution of detections in time

sfigure(1); %quiet figure that wont steal focus
clf
set(gcf,'color','w')%beige kills
histogram(data.counts_txy{1}(:,1),linspace(0,1,1e3))
xlabel('time(s)')
ylabel('counts')

%% Combining multiple shots
% now lets plot the same thing but for the first 10 shots combined
txy_multiple_shots=vertcat(data.counts_txy{1:10});

sfigure(1);
clf
set(gcf,'color','w')
[counts,edges] = histcounts(txy_multiple_shots(:,1),linspace(0,1,1e5));
edges=(edges(2:end)+edges(1:end-1))/2; %find the bin center
plot(edges,counts,'k')
xlabel('time(s)')
ylabel('counts')

%% FFT
%lets take the FFT of this histogram
sfigure(1)
set(gcf,'color','w')
fftout=fft_tx(edges,counts,'padding',3,'window','blackman');
semilogy(fftout(1,:),abs(fftout(2,:)))
xlabel('Frequency (Hz)')
ylabel('count modulation')

% That should give you a decent start
% for more advanced applications you will want to loop over each shot
