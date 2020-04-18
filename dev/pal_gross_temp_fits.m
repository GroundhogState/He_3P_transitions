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

%%
single_shot_plot = false;
shot_data=[];
shot_data.temp = nan(num_shots,2);
shot_data.temp_SE = nan(num_shots,2);
shot_data.pulse_data = cell(num_shots,1);
shot_data.ok_mask = zeros(num_shots,1);
shot_data.shot_num= ones(num_shots,1);
shot_data.num_counts= ones(num_shots,1);
shot_data.x_residual= ones(num_shots,1);
shot_data.y_residual= ones(num_shots,1);
x_refit_lvl = 0.2;
y_refit_lvl = 0.2;
x_spatial_min = 0.015;
y_spatial_min = 0.015;
x_spatial_max = 0.035;
y_spatial_max = 0.035;
min_shot_counts = 1e4;
% dead shots: 41, 48, 110
for i=1:num_shots
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

%             cli_header(2,'Fitting...');
            width_guess = mean(pulse_data.std);
            X = pulse_data.x_cens(1,:);
            x_mean_flux = mean(pulse_data.x_flux);
            x_flux_SE = std(pulse_data.x_flux)/sqrt(num_pulses);
            x_mask = abs(X) > x_spatial_min & abs(X) < x_spatial_max;

            Y = pulse_data.y_cens(1,:);
            y_mean_flux = mean(pulse_data.y_flux);
            y_flux_SE = std(pulse_data.y_flux)/sqrt(num_pulses);
            y_mask = abs(Y) > y_spatial_min & Y < y_spatial_max;

            gfun = @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2))+b(4); 
            coeff_names={'amp','mu','sig','offset'};
            x_guess=[mean(x_mean_flux),0,2*width_guess(2),0];
            x_mdlfit = fitnlm(X(x_mask),x_mean_flux(x_mask),gfun,x_guess);
            x_fit = x_mdlfit.Coefficients.Estimate;
            x_resid = x_mean_flux-gfun(x_fit,X);
            x_rel_err = x_resid./x_mean_flux;
            x_refit = abs(x_rel_err) < x_refit_lvl;
            x_mask =  x_refit;
            x_mdlfit = fitnlm(X(x_mask),x_mean_flux(x_mask),gfun,x_guess);
            x_fit = x_mdlfit.Coefficients.Estimate;
            x_SE = x_mdlfit.Coefficients.SE;x_resid = x_mean_flux-gfun(x_fit,X);
            x_rel_err = x_resid./x_mean_flux;
            x_temperature_val=(abs(x_fit(3))*opts.global.fall_velocity/opts.global.fall_time)^2 *const.mhe/const.kb;
            x_temperature_unc=x_temperature_val*2*x_SE(3)/abs(x_fit(3));
            
            
            y_guess=[mean(x_mean_flux),0,2*width_guess(2),0];
            y_mdlfit = fitnlm(Y(y_mask),y_mean_flux(y_mask),gfun,y_guess);
            y_fit = y_mdlfit.Coefficients.Estimate;
            y_resid = y_mean_flux-gfun(y_fit,Y);
            y_rel_err = y_resid./y_mean_flux;
            y_refit = abs(y_rel_err) < y_refit_lvl;
            y_mask =  y_refit;
            y_mdlfit = fitnlm(Y(y_mask),y_mean_flux(y_mask),gfun,y_guess);
            y_fit = y_mdlfit.Coefficients.Estimate;
            y_SE = x_mdlfit.Coefficients.SE;
            y_resid = y_mean_flux-gfun(y_fit,Y);
            y_rel_err = y_resid./y_mean_flux;            
            y_temperature_val=(abs(y_fit(3))*opts.global.fall_velocity/opts.global.fall_time)^2 *const.mhe/const.kb;
            y_temperature_unc=y_temperature_val*2*y_SE(3)/abs(y_fit(3));
%             cli_header(2,'Plotting...');

            %% Consider: Adding a bit here that adaptively redefines fit domain
            %% to include unfitted regions with small residuals?


            if single_shot_plot    
                stfig('Pulse data');
                clf
                subplot(2,2,1)
                hold on
                plot(pulse_data.x_cens(1,:),x_mean_flux,'k')
                plot(pulse_data.x_cens(1,:),x_mean_flux+x_flux_SE,'k:')
                plot(pulse_data.x_cens(1,:),x_mean_flux-x_flux_SE,'k:')
                plot(X(x_mask),x_mean_flux(x_mask),'r.')
                plot(pulse_data.x_cens(1,:),gfun(x_guess,pulse_data.x_cens(1,:)),'g:')
                plot(pulse_data.x_cens(1,:),gfun(x_fit,pulse_data.x_cens(1,:)),'g')
                ylim([10,max(2*x_mean_flux)])
                set(gca,'Yscale','log')
                title('Mean X profile')

                subplot(2,2,2)
                hold on
                plot(pulse_data.y_cens(1,:),y_mean_flux,'k')
                plot(pulse_data.y_cens(1,:),y_mean_flux+y_flux_SE,'k:')
                plot(pulse_data.y_cens(1,:),y_mean_flux-y_flux_SE,'k:')
                plot(Y(y_mask),y_mean_flux(y_mask),'r.')
                plot(pulse_data.y_cens(1,:),gfun(y_guess,pulse_data.y_cens(1,:)),'g:')
                plot(pulse_data.y_cens(1,:),gfun(y_fit,pulse_data.y_cens(1,:)),'g')
                ylim([10,max(2*y_mean_flux)])
                set(gca,'Yscale','log')
                title('Mean Y profile')

                subplot(2,2,3)
                plot(X,x_resid./x_mean_flux,'k.','MarkerSize',0.5)
                ylim([-1,1])
                title('X relative error')
                subplot(2,2,4)
                plot(Y,y_resid./y_mean_flux,'r.','MarkerSize',0.5)
                ylim([-1,1])
                title('Y relative erros')

                suptitle(sprintf('Single shot fit %u',i));
                drawnow()
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
    
    pulse_data.x_residual = x_mean_flux-gfun(x_fit,X);
    pulse_data.y_residual = y_mean_flux-gfun(y_fit,Y);
    shot_data.temp(i,:) = [x_temperature_val,y_temperature_val];
    shot_data.temp_unc(i,:) = [x_temperature_unc,y_temperature_unc];
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
plot(shot_data.shot_num(logical(shot_data.ok_mask)),X_temps,'k.')
plot(shot_data.shot_num(logical(shot_data.ok_mask)),Y_temps,'r.')
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
for ii=1:num_shots
    if shot_data.ok_mask(ii)
    px=plot(shot_data.pulse_data{ii}.x_cens(shot_data.x_mask),...
            shot_data.pulse_data{ii}.x_residual(shot_data.x_mask),'k.','MarkerSize',0.3);
    end
end
for ii=1:num_shots
    if shot_data.ok_mask(ii)
    py=plot(shot_data.pulse_data{ii}.x_cens(y_mask),...
            shot_data.pulse_data{ii}.y_residual(y_mask),'r.','MarkerSize',0.3);
    end
end
xlabel('X,Y (m)')
ylabel('Fit error')
legend([px,py],'X','Y')
title('Residuals');

subplot(2,2,4)
hold on
plot(shot_data.num_counts(shot_data.ok_mask),shot_data.temp(shot_data.ok_mask,1),'k.')
plot(shot_data.num_counts(shot_data.ok_mask),shot_data.temp(shot_data.ok_mask,2),'r.')
xlabel('Num counts')
ylabel('Temperature')
suptitle('Pulse-aligned thermometry')
