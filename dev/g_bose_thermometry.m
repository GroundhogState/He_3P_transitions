% 
%  Measured 5 microkelvin.
% % condensed frac = 1-(T/T_c)^3
% % So what's T_c?
% % Standard trap 
% omega = [425,425,45];
% 
% prop=bec_properties(omega,7.5e5,const.ahe_scat);
% prop.tc
% T=5e-6;
% 1-(T/prop.tc.non_interacting)^3
% 
% Thermal cloud expands like
% n_p(r,t) = n_p(r-p/m t). One finds
% % Pitaevskii eqn (10.40)
% n_T(r,t) = lambda_dB^-3 * g_(3/2)(z(r,t)) \prod 1/sqrt(1 + omega_k^2 t^2) 
% % where z(r,t) = exp(-beta V(r,t))
% % and V(r,t) = 1/2 m \sum omega_k^2/sqrt(1 + omega_k^2 t^2) r_k^2 plays the
% % role of the effective potentail, characteriing the shape of thermal
% % density
% % in the large-time limit this goes to
% % n_0,T (r,t) \rightarrow (m/t)^3 n_0,T(p)
% % eqn 10.31 - used in QD experiment - will have to try apply to this also,
% % darn

mcp_data = data.tdc;
mcp_data.all_ok = ones(size(mcp_data.N_atoms));
num_shots = length(mcp_data.N_atoms);
opts.global.fall_velocity = const.g0*opts.const.fall_time;
opts.global.fall_time = opts.const.fall_time;

opts.t0 = 0.4145; % Centre of first pulse (sec)
opts.pulsedt = .008; %time between pulses (sec)
opts.num_bins = [50,50,50];
opts.hist_2d = false;
opts.draw_plots = false;
opts.centre_bec = true;
opts.verbose = false;

min_flux = 1e2;

opts.xylim = [-0.03,.03;-.03,.03];
opts.pulse_twindow = .003; 

num_pulses = 150;

% %
% omega = [425,425,45];
% lambda_db = @(T) const.h/sqrt(2*pi*const.mhe*const.kb*T);
% omega_ho = 2*pi*geomean(omega);
% T = ((const.h/(A^(-1/3))/*const.mhe*omega_ho)^2)/(2*pi*const.mhe*const.kb)

g_tol = 1e-6;
n_p = @(p,v) p(1)*g_bose(exp(-(const.mhe*(v-p(2))).^2/(2*const.mhe*const.kb*p(3))),g_tol)+p(4);


A = 3e3;
T_guess=3e-7;
fit_guess = [A,0,T_guess,0];





% also have <p^2>_T = zeta(4)/zeta(3) * m k_B T/2 independent of direction

%%

single_shot_plot = true;

x_refit_lvl = 0.2;
y_refit_lvl = 0.2;
x_spatial_min = 0.015;
y_spatial_min = 0.015;
x_spatial_max = 0.035;
y_spatial_max = 0.035;
min_shot_counts = 1e4;


shot_data=[];
shot_data.temp = nan(num_shots,2);
shot_data.temp_SE = nan(num_shots,2);
shot_data.pulse_data = cell(num_shots,1);
shot_data.ok_mask = zeros(num_shots,1);
shot_data.shot_num= ones(num_shots,1);
shot_data.num_counts= ones(num_shots,1);
shot_data.x_residual= ones(num_shots,1);
shot_data.y_residual= ones(num_shots,1);
% dead shots: 41, 48, 110
for i=1:num_shots
    warning('off')
    if mod(i-1,10)==0
        cli_header('Working on shot %u...',i);
    end
%     try
        mask = zeros(size(data.tdc.N_atoms));
        mask(i) = 1;
        tshot = struct_mask(data.tdc,logical(mask));
        if tshot.N_atoms < min_shot_counts
            warning('Insufficient counts')
            shot_data.ok_mask(i) = 0;
        else
            pulse_data.x_cens = zeros(opts.num_bins(2));
            pulse_data.y_cens= zeros(opts.num_bins(3));
            pulse_data.x_flux= zeros(num_pulses,opts.num_bins(2));
            pulse_data.y_flux= zeros(num_pulses,opts.num_bins(3));
            pulse_data.txy = cell(num_pulses,1);
            pulse_data.num= zeros(num_pulses,1);
            pulse_data.cen= zeros(num_pulses,3);
            pulse_data.std= zeros(num_pulses,3);

            for pulse_num = 1:num_pulses
                opts.lims = [opts.t0+(pulse_num-1)*opts.pulsedt,opts.t0+(pulse_num)*opts.pulsedt;
                            opts.xylim(1,:);
                            opts.xylim(2,:)];
                h_data = show_txy_raw(tshot,opts);
                pulse_data.x_cens = h_data.centres{2};
                pulse_data.y_cens = h_data.centres{3};
                pulse_data.x_flux(pulse_num,:) = h_data.flux{2};
                pulse_data.y_flux(pulse_num,:) = h_data.flux{3};
                pulse_data.txy{pulse_num} = h_data.txy;
                pulse_data.num(pulse_num) = size(h_data.txy,1);
                pulse_data.cen(pulse_num,:) = h_data.pulse_cen;
                pulse_data.std(pulse_num,:) = h_data.pulse_std;
            end  % loop over pulses

            
            width_guess = mean(pulse_data.std);
            X = pulse_data.x_cens(1,:);
            x_mean_flux = mean(pulse_data.x_flux);
            x_flux_SE = std(pulse_data.x_flux)/sqrt(num_pulses);
    
            Y = pulse_data.y_cens(1,:);
            y_mean_flux = mean(pulse_data.y_flux);
            y_flux_SE = std(pulse_data.y_flux)/sqrt(num_pulses);

%             cli_header(2,'Fitting X...');
            v_x = X/opts.const.fall_time;
            x_mask = abs(v_x) > .01;
            x_mdl = fitnlm(v_x(x_mask),x_mean_flux(x_mask),n_p,fit_guess);
            x_coef = x_mdl.Coefficients.Estimate;
            x_SE = x_mdl.Coefficients.SE;

%             cli_header(2,'Fitting Y...');
            v_y = Y/opts.const.fall_time;
            y_mask = abs(v_y) > .035 & y_mean_flux > 100;
            y_mdl = fitnlm(v_y(y_mask),y_mean_flux(y_mask),n_p,fit_guess);
            y_coef = y_mdl.Coefficients.Estimate;
            y_SE = y_mdl.Coefficients.SE;

            %             cli_header(2,'Plotting...');

%             cli_header(1,'Temperatures fitted: %.3e, %.3e',x_coef(3),y_coef(3));
            if single_shot_plot    
                cli_header(2,'Plotting...');
                stfig('Pulse data');
                clf
                subplot(2,2,1)
                hold on
                plot(v_x,x_mean_flux,'k')
                plot(v_x(x_mask),x_mean_flux(x_mask),'b.')
                plot(v_x,n_p(fit_guess,abs(v_x)),'r:')
                plot(v_x,n_p(x_coef,abs(v_x)),'r')
                ylim([100,max(2*x_mean_flux)])
                xlabel('$v_x$')
                set(gca,'Yscale','log')
                title('Mean X profile')

                subplot(2,2,2)
                hold on
                plot(v_y,y_mean_flux,'k')
                plot(v_y(y_mask),y_mean_flux(y_mask),'b.')
%                 plot(v_y,n_p(x_coef,abs(v_y)),'r:')
                plot(v_y,n_p(y_coef,abs(v_y)),'r')
                set(gca,'Yscale','log')
                title('Mean Y profile')
                
                subplot(2,2,3)
                plot(v_x(x_mask),n_p(x_coef,abs(v_x(x_mask)))-x_mean_flux(x_mask),'k.')
                
                subplot(2,2,4)
                plot(v_y(y_mask),n_p(y_coef,abs(v_y(y_mask)))-y_mean_flux(y_mask),'r.')

                suptitle(sprintf('Single shot fit %u',i));
            end % single shot plot
            shot_data.ok_mask(i) = 1;
        end %if enough counts
        [msgstr, msgid] = lastwarn; %were there warnings in that attempt?
        if strcmp(msgid,'stats:nlinfit:IllConditionedJacobian')
            cli_header(2,'Omitting shot %u, poorly fitted\n',i);
            shot_data.ok_mask(i) = 0;
            warning('Restting warning...');
        end
        
%     catch
%         warning('Error in shot %u',i);
%         shot_data.ok_mask(i) = 0;
%     end %try shot
%     cli_header(3,'Done.');
    
    pulse_data.x_residual = x_mean_flux-n_p(x_coef,X);
    pulse_data.y_residual = y_mean_flux-n_p(y_coef,Y);
    shot_data.temp(i,:) = [x_coef(3),y_coef(3)];
    shot_data.temp_unc(i,:) = [x_SE(3),y_SE(3)];
    shot_data.amp(i,:) = [x_coef(1),y_coef(1)];
    shot_data.amp_unc(i,:) = [x_SE(1),y_SE(1)];
    shot_data.pulse_data{i} = pulse_data;
    shot_data.shot_num(i) = i;
    shot_data.x_mask = x_mask;
    shot_data.y_mask = y_mask;
    shot_data.num_counts(i) = mcp_data.N_atoms(i);
%     end
end %loop over shots
%%
shot_data.ok_mask = logical(shot_data.ok_mask);
idxs = logical(shot_data.ok_mask);
cli_header(1,'Plotting...');
stfig('PAL temp fits');
clf
subplot(2,2,1)
hold on
X_temps = 1e6*shot_data.temp(logical(shot_data.ok_mask),1);
Y_temps = 1e6*shot_data.temp(logical(shot_data.ok_mask),2);
plot(shot_data.shot_num(logical(shot_data.ok_mask)),X_temps,'k.','MarkerSize',3)
plot(shot_data.shot_num(logical(shot_data.ok_mask)),Y_temps,'r.','MarkerSize',3)
ylim([0,max(max([X_temps,Y_temps]))])
xlabel('Shot number')
ylabel('Temperature ($\mu$K)')
legend('X','Y')
title('Fitted temperatures ($\mu$K)')
subplot(2,2,2)
hold on
histogram((shot_data.temp(:,1)-shot_data.temp(:,2))./nanmean(shot_data.temp,2),30,'FaceColor',[.1,.1,.7],'FaceAlpha',0.4)
xlabel('Relative err (x-y) / 2(x+y)')
ylabel('Num shots')
title('Relative difference in X,Y temperatures')



subplot(2,2,3)
hold on
plot(shot_data.num_counts(shot_data.ok_mask),shot_data.temp(shot_data.ok_mask,1),'k.','MarkerSize',0.5)
plot(shot_data.num_counts(shot_data.ok_mask),shot_data.temp(shot_data.ok_mask,2),'r.','MarkerSize',0.5)
xlabel('Num counts')
ylabel('Temperature')
suptitle('Pulse-aligned thermometry')


