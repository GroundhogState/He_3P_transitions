c = init_constants;
data_stash = "C:\Users\jaker\GoogleDrive\HEBEC\Projects\core_Projects\He_3P_transitions\dev";
fname = "ZeemanShifts_inMHz";
% Two-peak results STAGE 1:
% All peaks: [-12.87,-17.15,-47.90,14.65]
% relative to 744396208.36
% Two-peak results STAGE 2:
% All peaks: [-0.07,-3.69,-25.99,16.28]
% relative to 744396208.36


f_23P2_53D1 = 744396511;
f_23P2_53D2 = 744396227;
f_23P2_53D3 = 744396208;
f_dref= 744396208.36;
f_pred = [f_23P2_53D1,f_23P2_53D2,f_23P2_53D3];
E3 = f_23P2_53D1 - f_23P2_53D3;
E2 = f_23P2_53D2 - f_23P2_53D3;
f_23P = @(B) B*c.mu*2*lande_g(1,2,2)/(1e6*c.h);
cell_shift = -1.9; %MHz
% expdata{1}.peaks = [744396158.56 744396191.07 744396221.12,744396451.51]; % stage 1
% expdata{2}.peaks = [744396180.504 744396204.973 744396222.922,744396475.92]; % stage 2
expdata.peaks = [f_dref+[-47.90,-17.15,-12.87,14.65],744396451.51; % stage 1
                f_dref+[-25.99,-3.69,-0.07,16.28],744396475.92]; % stage 2
expdata.field= [18.25 11.43];

% Energy levels of the 5 ^3D_j system with mJ=1 and #3 of the set. 
% First row is 0 followed by list of j=3->j=2 splittings (E2). 
% First column is 0 followed by magnetic fields in G.
% Set up the hamiltonian
% %

E_1 = 0; %MHz
E_2 = -280; %MHz
E_3 = -300; %MHz
E_lvl = [E_1,E_2,E_3];
% H_init = Hamiltonian(E_lvl,B,c);
num_B = 20;
plot_fields = linspace(0.001,20 ,num_B);
evals = zeros(num_B,5);
% evals_opt = zeros(num_B,5);
for i = 1:num_B
    evals(i,:) = sort(eig(Hamiltonian(E_lvl,plot_fields(i),c)));
end
% Optimize the levels
opt_fun = @(x) split_fun(x,expdata,c);
opt_lvls = fminsearch(opt_fun,E_lvl);
evals_opt = zeros(num_B,5);
for i = 1:num_B
    evals_opt(i,:) = eig(Hamiltonian(opt_lvls,plot_fields(i),c));
end
% %
% Ok, there is a way to allocate these curves: Start with the initial
% eigenspectrum. Then diagonalize at the next field value, and for each
% eigenvector(t), find the eigenvector(t+1) that maximizes the dot product.
% That is, take the eigenvector matrices and multiply them!
% Then mask this to a permutation matrix.
fprintf('Drawing fields...\n')
all_evals = zeros(num_B,5);
all_evecs = zeros(num_B,5,5);
% all_proj =  zeros(num_B,5,5);
eval_labels = zeros(num_B,5);



%% Make plots
ev_or_spec = 1;
data_colours = [226,90,42;
                45,225,245]/255;
f1=stfig('Constraining interval with data');
clf;
hold on
Y1 = expdata.peaks(1,:)-f_23P2_53D1;
Y2 = expdata.peaks(2,:)-f_23P2_53D1;
X1 = expdata.field(1)*ones(size(Y1));
X2 = expdata.field(2)*ones(size(Y1));
Yerr = 2.2*ones(size(Y1));
Xerr = 0.3*ones(size(Y1));
% th2=plot(plot_fields,(evals_opt(:,2))/1e6-f_23P(plot_fields)','LineWidth',1.5);

th1=plot(plot_fields,(evals_opt(:,2))/1e6-ev_or_spec*f_23P(plot_fields)',':','LineWidth',1.5,'Color',[0.5,0.5,0.5]);
th2=plot(plot_fields,(evals_opt(:,4))/1e6-ev_or_spec*f_23P(plot_fields)','-','LineWidth',1.5,'Color',[0.5,0.5,0.5]);
th3=plot(plot_fields,(evals_opt(:,3))/1e6-ev_or_spec*f_23P(plot_fields)','--','LineWidth',1.5,'Color',[0.5,0.5,0.5]);
th4=plot(plot_fields,(evals_opt(:,5))/1e6-ev_or_spec*f_23P(plot_fields)','-.','LineWidth',1.5,'Color',[0.5,0.5,0.5]);
e1=errorbar(X1,Y1,Yerr,Yerr,Xerr,Xerr,'.','MarkerSize',25,'MarkerFaceColor',data_colours(1,:));
e2=errorbar(X2,Y2,Yerr,Yerr,Xerr,Xerr,'.','MarkerSize',25,'MarkerFaceColor',data_colours(2,:));
legend([th1,th2,th3,th4],{'$|21\rangle$','$|22\rangle$','$|31\rangle$','$|32\rangle$'},'Location','SouthWest')
if ev_or_spec==0
  ylim([-370,-150])  
else
    ylim([-370,-270])
end
set(gca,'FontSize',20)
title('Fitting $5^{3\!}D$ levels','FontSize',35)
xlabel('Magnetic field strength (G)','FontSize',30)
ylabel(sprintf('$\\nu$ - %u (MHZ)',f_23P2_53D1),'FontSize',30)

% 
% cli_header(0,'Results of optimization:');
% for i=1:3
%     cli_header(1,'J=%u level %.2f',i,f_23P2_53D1+opt_lvls(i));
% end

% cli_header(2, 'Lines at:');
% sort(eig(Hamiltonian(opt_lvls,expdata.field(1),c)))/1e6'-0*f_23P(expdata.field(1))
% sort(eig(Hamiltonian(opt_lvls,expdata.field(2),c)))/1e6'-0*f_23P(expdata.field(2))
% cli_header(1,'Relative to %.2f',f_23P2_53D1);

cli_header(0,'All done.');
% %
savedir = 'C:\Users\jaker\GoogleDrive\HEBEC\Thesis\fig\Spectroscopy';
imname = sprintf('solving_5D_splitting');
filename1 = fullfile(savedir,imname);
            saveas(f1,[filename1,'.fig']);
            saveas(f1,[filename1,'.png']);
            saveas(f1,[filename1,'.eps']);
            saveas(f1,[filename1,'.svg']);
cli_header(1,'Done.');



%%

function H = Hamiltonian(E_lvl,B,c)
    U = [sqrt(1/10) ,-sqrt(3/10)	,   0		,  sqrt(3/5)	,    0	;         %\J mJ = 11
        -sqrt(1/2)	, sqrt(1/6)	,	0 		,  sqrt(1/3)	,    0	;             %\J mJ = 21
        0			,	0			,-sqrt(1/3),  	0  			,   sqrt(2/3) 	; %\J mJ = 22
        sqrt(2/5)  ,	sqrt(8/15)	,	 0		,  sqrt(1/15)  ,    0 	;         %\J mJ = 31
        0			,	0			,sqrt(2/3) ,  0   			,    sqrt(1/3)];  %\J mJ = 32
    H_alpha = diag([2,1,3,0,2]);
    H_0 = diag([E_lvl(1),E_lvl(2),E_lvl(2),E_lvl(3),E_lvl(3)])*1e6;
    H_v = c.mu/c.h*U'*H_alpha*U;
    H = H_0 + B*H_v;
end
    

function sqr_err = split_fun(E_lvl, expdata,c)
    f_23P2_53D1 = 744396511;
    H_1 = Hamiltonian(E_lvl,expdata.field(1),c);
    H_2 = Hamiltonian(E_lvl,expdata.field(2),c);
    f_23P = @(B) B*c.mu*2*lande_g(1,2,2)/(1e6*c.h);
    freqs_stage = zeros(2,length(H_1));
    exp_splittings = zeros(2,length(expdata.peaks(1,:)));
    freqs_stage(1,:) = sort(eig(H_1)/1e6 - f_23P(expdata.field(1)));
    freqs_stage(2,:) = sort(eig(H_2)/1e6 - f_23P(expdata.field(2)));
    exp_splittings(1,:) = sort(expdata.peaks(1,:)-f_23P2_53D1);
    exp_splittings(2,:) = sort(expdata.peaks(2,:)-f_23P2_53D1);
    split_diff = exp_splittings - freqs_stage;
    sqr_err = sum(sum((split_diff).^2));
end


%%

%%
% function square_err = errfunc(E_vec,fields,measured_freqs)
%     predictions = 
% end
% 
% function squerr = err_vs_splitting(field_idxs,split_idx,levels,expdata)
%     d_theory_1 = squeeze(levels(2,2,field_idxs(1),split_idx))-squeeze(levels(1,1,field_idxs(1),split_idx));
%     d_theory_2 = squeeze(levels(2,2,field_idxs(2),split_idx))-squeeze(levels(1,1,field_idxs(2),split_idx));
%     d_exp_1 = expdata{1}.peaks(3)-expdata{1}.peaks(1); % high field!
%     d_exp_2 = expdata{2}.peaks(3)-expdata{2}.peaks(1); % low field!
%     % Plot the overlapping peaks
%     % plot(fields, squeeze(levels(1,2,:,split_idx))-f_lower)
%     % plot(fields, squeeze(levels(2,1,:,split_idx))-f_lower)
%     % Well we can do this the dumb way...
%     squerr = (d_theory_2-d_exp_1).^2+(d_theory_1-d_exp_2).^2;
% end

%% Original work with Danny's tables

% ax_ref = csvread(fullfile(data_stash,"ZeemanShifts_inMHz_mJ_1_ind_1.csv"),2);
% splittings = ax_ref(1,2:end);
% fields = ax_ref(2:end,1);
% levels = zeros(3,3,101,101);
% for mJ = 1:3
%     for idx = 1:4-mJ
%         levels(mJ,idx,:,:) = csvread(fullfile(data_stash,sprintf("ZeemanShifts_inMHz_mJ_%u_ind_%u.csv",mJ,idx)),3,1);
%     end
% end
% 
% % Vary this and find min distance to peaks
% split_idx= find(splittings == 21);
% field_idx_1 = find(fields == 11);
% field_idx_2 = find(fields == 18);
% % Plot the difference from the f_23P2_53D2 peak at 744396227 MHz
% % Adding the theory bits... So we have the splittings relative to J=3 level
% % Ok, well then just change the f_23P2_53D1, right? Hm, still a sign problem.
% f_23P2_53D1 = f_23P2_53D3;
% f_lower = fields*c.mu*2*lande_g(1,2,2)/(1e6*c.h);
% f_idxs = [field_idx_1,field_idx_2];
% stfig('Splitting calculations');
% clf;
% hold on
% squerrs = err_vs_splitting(f_idxs,1:100,levels,expdata);
% 
% subplot(2,2,[1,2])
% hold on
% % plot empirical differences
% plot(expdata{1}.field,expdata{1}.peaks(3)-f_23P2_53D1,'kx')
% plot(expdata{1}.field,expdata{1}.peaks(1)-f_23P2_53D1,'kx')
% plot(expdata{2}.field,expdata{2}.peaks(3)-f_23P2_53D1,'kx')
% plot(expdata{2}.field,expdata{2}.peaks(1)-f_23P2_53D1,'kx')
% % And the overlapping peaks
% plot(expdata{2}.field,expdata{2}.peaks(2)-f_23P2_53D1,'rx')
% plot(expdata{1}.field,expdata{1}.peaks(2)-f_23P2_53D1,'rx')
% plot(expdata{1}.field,D1_data(1)-f_23P2_53D1,'bx')
% plot(expdata{2}.field,D1_data(2)-f_23P2_53D1,'bx')
% plot(expdata{1}.field,744396236.68-f_23P2_53D1,'ko')
% cols = {'k','b','r'};
% for mJ = 1:3
%     for idx = 1:4-mJ
%         plot(fields,squeeze(levels(mJ,idx,:,split_idx))-f_lower-12,cols{mJ})
%     end
% end
% ylim([-100,120])
% ylabel('df')
% xlabel('Field strength (G)')
% title('Predicted lines for a field-free splitting of 21MHz')
% 
% subplot(2,2,[3,4])
% hold on
% plot(splittings(1:60),squerrs(1:60),'k')
% xlabel('$D_2-D_3$ splitting')
% ylabel('Peak difference error$^2$')
% title('Squared error Theory-Expt')








