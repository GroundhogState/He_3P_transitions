%Analysis of 5^3D states
this_folder = fileparts(which(mfilename));
core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
addpath(genpath(this_folder));
addpath(genpath(core_folder));
% Physical constants
const = init_constants();
% Table of field-free freqs defined here for readability
% Read N_(2S+1)L_J 
% The 412nm benchmark
const.f_table.g_2_3P_2.e_5_3S_1 = 727.3032446e12;
% Misc transitions - what do the stars mean?
const.f_table.g_2_3P_2.e_5_3P_0 = 1e9*const.c/404.628937550957;
const.f_table.g_2_3P_2.e_5_3P_1 = 1e9*const.c/404.629844755577;
const.f_table.g_2_3P_2.e_5_3P_2 = 1e9*const.c/404.629918705477; 
% Historically controversial transitions
    % Updated c.f. Drake's email
const.f_table.g_2_3P_2.e_5_3D_3 = 744396208.36e6;%744.39620968e12;
const.f_table.g_2_3P_2.e_5_3D_2 = 744396227.58e6;% 744.39622889e12;
const.f_table.g_2_3P_2.e_5_3D_1 = 744396511.14e6 ;%744.39651246e12; 
% Singlet-triplet transitions
const.f_table.g_2_3P_2.e_5_1S_0 = 1e9*const.c/406.8886971706; % Can't observe from our pump state :( 
const.f_table.g_2_3P_2.e_5_1P_1 = 1e9*const.c/402.322271224483;  % Should be visible with sigma-
const.f_table.g_2_3P_2.e_5_1D_2 = 744430343.14e6;% Spotted; should be able to get two lines


%Fitted values for the 5^3D's
const.f_table.g_2_3P_2.e_5_3D_3 = 744.39620836e12;
const.f_table.g_2_3P_2.e_5_3D_2 = 744.39622758e12;
const.f_table.g_2_3P_2.e_5_3D_1 = 744.39651114e12;
g_state = '2_3P_2_2';
% ITC 1         ITC 2

% data ={744396158380000,744396180286000,'5_3D_3_1'
%            744396191764000,744396205362000,'5_3D_2_1'
%            744396220906000,744396225174000,'5_3D_2_2'
%            744396235218000,744396225174000,'5_3D_3_3'
%https://www.nist.gov/sites/default/files/documents/srd/jpcrd382009565p.pdf
data ={744396448622000,'ITC_1','5_3D_1_1'
       744396473900000,'ITC_2','5_3D_1_1'
       727303233450000,'ITC_2','5_3S_1_1'};
transitions = {};
B=[17.86,10.5];%the three stages
% figure(404);
% clf
B_plot = B.*ones(size(data,1),2);
for ii =1:size(data,1)
    e_state = data{ii,3};
    trans_name = ['t',e_state(1:6)];
    if ~isfield(transitions,(trans_name))
        transitions.(trans_name) = [];
    end
%     hold on
%     plot(B_plot(ii,:),cell2mat(data(ii,1:2)),'x-')
    if strcmp(data{ii,2},'ITC_1')
        B = 17.86;
    else
        B = 10.5;
    end

    
    f_measure = data{ii,1};
    lines = zeeman_state2state(B,g_state,e_state,const);
    f_adj = f_measure-lines.df;
    transitions.(trans_name) = [transitions.(trans_name),f_adj];
%     f_0 = nanmean(f_adj);
%     f_mdl = fitlm(B,f_measure);
end
% xlabel('Magnetic field strength (G)')
% ylabel('Frequency (Hz)')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% WRITE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwtext(' RESULTS ')
header({1,'Zeeman Corrected Transitions'})
fields = fieldnames(transitions);
for idx = 1:length(fields)
   current_trans = transitions.(fields{idx});
    fprintf('Transition Frequency %s             %.2f MHz\n',fields{idx},nanmean(current_trans)./1e6)
end
% fprintf('Peak width                       %.3f±(%.3f) MHz\n',data.gaus_fit(3),data.gaus_fit(4))
% fprintf('Transition Wavelength            %.6f±(%.7f) nm\n',data.gaus_fit(5),data.gaus_fit(6))
% header({1,'Lorentzian fit'})
% fprintf('Transition Frequency             %.2f±(%.3f) MHz\n',data.fits.lorz_fit(1),data.lorz_fit(2))
% fprintf('Peak width                       %.3f±(%.3f) MHz\n',data.lorz_fit(3),data.lorz_fit(4))
% fprintf('Transition Wavelength            %.6f±(%.7f) nm\n',data.lorz_fit(5),data.lorz_fit(6))
%%
function const = init_constants()

        const.mu = 9.27e-28; %J/G
        const.h = 6.63e-34;
        const.hbar = const.h/(2*pi);
        const.f_mu = const.mu/const.h;
        const.w_mu = const.mu/const.hbar;
        const.c = 299792458;
        % Notation & lookup
        const.terms = {'S','P','D','F','G'};

end