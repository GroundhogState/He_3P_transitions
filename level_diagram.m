%% Making a level diagram
% Set up
% clear all
profile on
% this_folder = fileparts(which(mfilename));
% addpath(genpath(this_folder));
% core_folder = fullfile(fileparts(this_folder),'Core_BEC_Analysis\');
% addpath(genpath(core_folder));
% zeeman_folder = fullfile(fileparts(this_folder),'Zeeman_splitting\');
% addpath(genpath(zeeman_folder));

opts = master_transition_config('null');
const = opts.const;
% Define starting information
cooling_wl = 1083.331e-9; %nm
cool_gap = const.c/cooling_wl; %Hz
metastable_gap = 19.8*const.q/const.h; %Hz

e_states = {'5^3S_1','5^3D_1','5^3D_2','5^3D_3','5^1D_2'};
% e_states = {'5^3S_1','5^3D_1','5^1D_2'};
labels = {'1^1S_0','2^3S_1','2^3P_2',e_states{:}};
lvls.labels = cellfun(@(x) x(1:4), labels,'uni',0);
lvls.states = cellfun(@(x) strrep(x,'^','_'),labels,'uni',0);

transitions.ends{1} = {'2_3S_1','2_3P_2'};
transitions.ends{2} = {'2_3P_2','5_3S_1'};
transitions.ends{3} = {'2_3P_2','5_3D_1'};
transitions.ends{4} = {'2_3P_2','5_3D_2'};
transitions.ends{5} = {'2_3P_2','5_3D_3'};
transitions.ends{6} = {'2_3P_2','5_1D_2'};

% Presentation options
l_spacing = 1;
l_width = 0.8;
l_jitter = 0.75;
lbl_offset = 0.3;
lv_labelsize = 12;
sector_labelsize = 18;


%% MAKE FIGURE
ntrans = numel(transitions.ends);
% Calculate the levels
for nstate = 1:numel(e_states)
    s_upper = lvls.states{nstate+3};
    Ls(nstate) = find(strcmp(s_upper(4),const.terms))-1;
    gaps(nstate) = const.f_table.('g_2_3P_2').(['e_',s_upper])+cool_gap+metastable_gap; %Hz
    Ss(nstate) = 0.5*(str2num(s_upper(3))-1);
end

lvls.levels = [0,metastable_gap,cool_gap+metastable_gap,gaps]*const.h/const.q; %eV
lvls.Ls = [0,0,1,Ls];
lvls.S = [0,1,1,Ss];
Lmax = max(lvls.Ls);

nstates = numel(lvls.states);

figure(29746);
clf;
set(gcf,'color','w');
sector_split = ((Lmax+0.5)*l_spacing);
plot(sector_split*[1,1],[0,25],'k:') % divide the spin sectors
hold on
% Plot the levels
for n = 1:nstates
    Y = lvls.levels(n);
    L_low = lvls.Ls(n);
    spin_shift_low = lvls.S(n)*(1+Lmax)*l_spacing;
    x0 = L_low*l_spacing-l_width/2+spin_shift_low;
    x1 = L_low*l_spacing+l_width/2+spin_shift_low;
    plot([x0,x1],(Y)*[1,1],'k')
    lvls.X(n,:) = [x0,x1];
    text(x0+0.3,Y+lbl_offset,lvls.labels{n},'FontSize',lv_labelsize)
end

% Plot the transition lines
for m=1:ntrans
   s_low = find(strcmp(transitions.ends{m}{1},lvls.states));
   L_low = lvls.Ls(s_low);
   spin_shift_low = lvls.S(s_low)*(1+Lmax)*l_spacing;
   x0 = L_low*l_spacing+spin_shift_low;
   
   s_hi = find(strcmp(transitions.ends{m}{2},lvls.states));
   L_hi = lvls.Ls(s_hi);
   spin_shift_hi= lvls.S(s_hi)*(1+Lmax)*l_spacing;
   x1 = L_hi*l_spacing+spin_shift_hi;
   
   if x0>x1 % moving left
       x0 = x0 - 0.5*l_jitter*l_width;
       x1 = x1 + 0.5*l_jitter*l_width;
   else
       x0 = x0 + 0.5*l_jitter*l_width;
       x1 = x1 - 0.5*l_jitter*l_width;
   end
   
   y0 = lvls.levels(s_low);
   y1 = lvls.levels(s_hi);
   
   plot([x0,x1],[y0,y1],'k')
    
end
%%
text(0.3*l_spacing,18.5,'Orthohelium \uparrow\downarrow','FontSize',sector_labelsize)
text(1*l_spacing+sector_split+1,18.5,'Parahelium \uparrow\uparrow','FontSize',sector_labelsize)


annotation('textarrow',[0.2,0.2],[0.3,0.15],'String','Ground state 20eV that way!','FontSize',12)
text([4,4],[20,20],'Cooling transition 1083.331nm')
xlim([-1,7])
ylim([-1,26])

suptitle('Helium level diagram')

box off
profile off
% profile viewer


%% 5^3D sublevels


slv = lvls.levels(5:7);
ymax = max(slv);
figure(245)
clf
set(gcf,'color','w');
for i=1:3
    Y = -1e6*(ymax-slv(i));
    x0 = L_low*l_spacing-l_width/2+spin_shift_low;
    x1 = L_low*l_spacing+l_width/2+spin_shift_low;
    plot([x0,x1],(Y)*[1,1],'k')
    hold on
end
box off
yticks(-1.2:.3:0)
ylim([-1.4,.1])
ylabel('\Delta E (\mu eV)')
xticks([])


%%
P_wnum = 100*[169087.8309,169086.8430,169086.7666]; %m^-1
P_lvl = const.hbar*const.c*P_wnum/const.q;
figure(246)
clf
set(gcf,'color','w');
for i=1:3
    Y = -1e6*(max(P_lvl)-P_lvl(i));
    x0 = L_low*l_spacing-l_width/2+spin_shift_low;
    x1 = L_low*l_spacing+l_width/2+spin_shift_low;
    plot([x0,x1],(Y)*[1,1],'k')
    hold on
end
box off
yticks(-20:10:0)
% ylim([-1.4,.1])
ylabel('\Delta E (\mu eV)')
xticks([])



