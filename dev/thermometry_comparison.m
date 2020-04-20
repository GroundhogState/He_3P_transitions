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
function shot_data = thermometry_comparison(data,opts)
    mcp_data = data.tdc;
    mcp_data.all_ok = ones(size(mcp_data.N_atoms));
    num_shots = length(mcp_data.N_atoms);
    opts.global.fall_velocity = 9.796*opts.const.fall_time;
    opts.global.fall_time = opts.const.fall_time;

    opts.t0 = 0.4145; % Centre of first pulse (sec)
    opts.pulsedt = .008; %time between pulses (sec)
    opts.num_bins = [50,50,50];
    opts.hist_2d = false;
    opts.draw_plots = false;
    opts.centre_bec = true;
    opts.verbose = false;



    opts.xylim = [-0.03,.03;-.03,.03];
    opts.pulse_twindow = .003; 



    % %
    % omega = [425,425,45];
    % lambda_db = @(T) const.h/sqrt(2*pi*const.mhe*const.kb*T);
    % omega_ho = 2*pi*geomean(omega);
    % T = ((const.h/(A^(-1/3))/*const.mhe*omega_ho)^2)/(2*pi*const.mhe*const.kb)

    g_tol = 1e-5;
    n_p = @(p,v) p(1)*g_bose(exp(-(const.mhe*(v-p(2))).^2/(2*const.mhe*const.kb*p(3))),g_tol)+p(4);
    gfun = @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3)))+b(4); 

    A = 3e3;
    T_guess=3e-7;
    fit_guess = [A,0,T_guess,0];





    % also have <p^2>_T = zeta(4)/zeta(3) * m k_B T/2 independent of direction

    %%
    g0 = 9.796;
    mhe = 6.6433e-27;
    kb=1.380e-23;
    % show_txy_raw
    min_flux = 1e4;
    num_pulses = 150;
    single_shot_plot = false;
    do_refits = true;
    min_shot_counts = 1e2;

    x_refit_lvl = 0.15;
    y_refit_lvl = 0.15;
    x_spatial_min = 0.015;
    y_spatial_min = 0.015;
    x_spatial_max = 0.035;
    y_spatial_max = 0.035;

    fit_options = statset('TolFun',1e-6);

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
        try
            mask = zeros(size(data.tdc.N_atoms));
            mask(i) = 1;
            tshot = struct_mask(data.tdc,logical(mask));
            if tshot.N_atoms < min_shot_counts
                warning('on')
                warning('Insufficient counts')
                warning('off')
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
                v_x = X/opts.const.fall_time;

                Y = pulse_data.y_cens(1,:);
                y_mean_flux = mean(pulse_data.y_flux);
                y_flux_SE = std(pulse_data.y_flux)/sqrt(num_pulses);
                v_y = Y/opts.const.fall_time; 

                x_mask = abs(v_x) > .01;
                x_mdl = fitnlm(v_x(x_mask),x_mean_flux(x_mask),n_p,fit_guess,...
                            'Options',fit_options);
                x_coef = x_mdl.Coefficients.Estimate;
                x_SE = x_mdl.Coefficients.SE;
                x_poly_rel_err = (n_p(x_coef,v_x)-x_mean_flux)./x_mean_flux;

                 % % Fit Gaussian functions
                x_guess_gauss=[mean(x_mean_flux),0,.0007,0];
                x_mdlfit_gauss = fitnlm(v_x(x_mask),x_mean_flux(x_mask),gfun,x_guess_gauss,...
                            'Options',fit_options);
                x_coef_gauss = x_mdlfit_gauss.Coefficients.Estimate;
                x_resid_gauss = x_mean_flux-gfun(x_coef_gauss,X);
                x_gauss_rel_err = (gfun(x_coef_gauss,v_x)-x_mean_flux)./x_mean_flux;   

                %re-fit
                x_refit_mask = abs(x_gauss_rel_err) < x_refit_lvl | abs(x_poly_rel_err) < x_refit_lvl;
                if sum(x_refit_mask) > sum(x_mask) && do_refits;
                    % New things to fit
                    x_mask = x_refit_mask & abs(v_x) > .01;
                    x_mdl = fitnlm(v_x(x_mask),x_mean_flux(x_mask),n_p,fit_guess,...
                                'Options',fit_options);
                    x_coef = x_mdl.Coefficients.Estimate;
                    x_SE = x_mdl.Coefficients.SE;
                    x_poly_rel_err = (n_p(x_coef,v_x)-x_mean_flux)./x_mean_flux;

                     % % Fit Gaussian functions
                    x_guess_gauss=[mean(x_mean_flux),0,.0007,0];
                    x_mdlfit_gauss = fitnlm(v_x(x_mask),x_mean_flux(x_mask),gfun,x_guess_gauss,...
                                'Options',fit_options);
                    x_coef_gauss = x_mdlfit_gauss.Coefficients.Estimate;
                    x_resid_gauss = x_mean_flux-gfun(x_coef_gauss,X);
                    x_gauss_rel_err = (gfun(x_coef_gauss,v_x)-x_mean_flux)./x_mean_flux;
                end            

                % Fitting Y
                % Polylog
                y_mask = abs(v_y) > 0.034 & y_mean_flux > 100;
                y_mdl = fitnlm(v_y(y_mask),y_mean_flux(y_mask),n_p,fit_guess,...
                        'Options',fit_options);
                y_coef = y_mdl.Coefficients.Estimate;
                y_SE = y_mdl.Coefficients.SE;
                y_poly_rel_err = (n_p(y_coef,v_y)-y_mean_flux)./y_mean_flux;       
                % Gaussian
                y_guess_guess=[mean(x_mean_flux),0,.0007,0];
                y_mdlfit_gauss = fitnlm(v_y(y_mask),y_mean_flux(y_mask),gfun,y_guess_guess,...
                        'Options',fit_options);
                y_coef_gauss = y_mdlfit_gauss.Coefficients.Estimate;
                y_resid_gauss = y_mean_flux-gfun(y_coef_gauss,Y);
                y_gauss_rel_err = (gfun(y_coef_gauss,v_y)-y_mean_flux)./y_mean_flux;

                %re-fit
                y_refit_mask = (abs(y_gauss_rel_err) < y_refit_lvl | abs(y_poly_rel_err) < y_refit_lvl) & abs(v_y) > 0.01;
                if sum(y_refit_mask) > sum(y_mask) && do_refits
                    % New things to fit
                    y_mask = y_refit_mask;
                    y_mdl = fitnlm(v_y(y_mask),y_mean_flux(y_mask),n_p,fit_guess,...
                                'Options',fit_options);
                    y_coef = y_mdl.Coefficients.Estimate;
                    y_SE = y_mdl.Coefficients.SE;
                    y_poly_rel_err = (n_p(y_coef,v_y)-y_mean_flux)./y_mean_flux;       
                    % Gaussian
                    y_guess_guess=[mean(x_mean_flux),0,.0007,0];
                    y_mdlfit_gauss = fitnlm(v_y(y_mask),y_mean_flux(y_mask),gfun,y_guess_guess,...
                                'Options',fit_options);
                    y_coef_gauss = y_mdlfit_gauss.Coefficients.Estimate;
                    y_resid_gauss = y_mean_flux-gfun(y_coef_gauss,Y);
                    y_gauss_rel_err = (gfun(y_coef_gauss,v_y)-y_mean_flux)./y_mean_flux;
                end


                polylog_r_squared = [x_mdl.Rsquared.Ordinary,y_mdl.Rsquared.Ordinary];
                gauss_r_squared = [x_mdlfit_gauss.Rsquared.Ordinary,y_mdlfit_gauss.Rsquared.Ordinary];
                x_temp_gauss =(abs(x_coef_gauss(3))) *const.mhe/const.kb;
                y_temp_gauss = (abs(y_coef_gauss(3))) *const.mhe/const.kb;

    %             [~, msgid] = lastwarn; %were there warnings in that attempt?
    %             if strcmp(msgid,'stats:nlinfit:IllConditionedJacobian')
    %                 cli_header(2,'Shot %u Gaussian poorly fitted\n',i);
    %     %             shot_data.ok_mask(i) = 0;
    %                 warning('Restting warning...');
    %             end


                pulse_data.x_residual = x_mean_flux-n_p(x_coef,X);
                pulse_data.y_residual = y_mean_flux-n_p(y_coef,Y);
                shot_data.temp(i,:) = [x_coef(3),y_coef(3)];
                shot_data.temp_unc(i,:) = [x_SE(3),y_SE(3)];
                shot_data.polylog_r_squared(i,:) = polylog_r_squared;
                shot_data.gauss_r_squared(i,:) = gauss_r_squared;
                shot_data.gauss_temp(i,:) = [x_temp_gauss,y_temp_gauss];
                shot_data.gauss_temp_unc(i,:) = [x_SE(3),y_SE(3)];
                shot_data.amp(i,:) = [x_coef(1),y_coef(1)];
                shot_data.amp_gauss(i,:) = [x_coef_gauss(1),y_coef_gauss(1)];
                shot_data.amp_unc(i,:) = [x_SE(1),y_SE(1)];
                shot_data.pulse_data{i} = pulse_data;
                shot_data.shot_num(i) = i;
                shot_data.x_mask = x_mask;
                shot_data.y_mask = y_mask;
                shot_data.num_counts(i) = mcp_data.N_atoms(i);

                if single_shot_plot   
                    cli_header(2,'Plotting...');
                    stfig('Pulse data');
                    clf
                    subplot(2,2,1)
                    hold on
                    plot(v_x,x_mean_flux,'k:')
                    plot(v_x(x_mask),x_mean_flux(x_mask),'k*')
                    plot(v_x(x_refit_mask),x_mean_flux(x_refit_mask),'ko')
                    plot(v_x,n_p(x_coef,abs(v_x)),'r')
                    plot(v_x,gfun(x_coef_gauss,v_x),'b')
                    ylim([100,max(2*x_mean_flux)])
                    legend('Flux','Fit domain','Polylog','Gaussian')
                    xlabel('$v_x$')
                    ylabel('Flux')
                    set(gca,'Yscale','log')
                    title('Mean X profile')

                    subplot(2,2,2)
                    hold on
                    plot(v_y,y_mean_flux,'k:')
                    plot(v_y(y_mask),y_mean_flux(y_mask),'k*')
                    plot(v_y(y_refit_mask),y_mean_flux(y_refit_mask),'ko')
                    plot(v_y,n_p(y_coef,abs(v_y)),'r')
                    plot(v_y,gfun(y_coef_gauss,v_y),'b')
                    xlabel('$v_y$')
                    ylabel('Flux')
                    set(gca,'Yscale','log')
                    title('Mean Y profile')

                    subplot(2,2,3)
                    hold on
    %                 plot(v_x(x_mask),(n_p(x_coef,abs(v_x(x_mask)))-x_mean_flux(x_mask))./x_mean_flux(x_mask),'ro')
    %                 plot(v_x(x_mask),(gfun(x_coef_gauss,(v_x(x_mask)))-x_mean_flux(x_mask))./x_mean_flux(x_mask),'bx')

                    plot(v_x,(n_p(x_coef,abs(v_x))-x_mean_flux)./x_mean_flux,'ro')
                    plot(v_x,(gfun(x_coef_gauss,(v_x))-x_mean_flux)./x_mean_flux,'bx')
                    legend('Polylog','Gaussian','Location','NorthEast')
                    ylim([-1,1])
                    xlabel('$v_x$')
                    ylabel('$\Delta$')
                    title('$v_x$ residuals')

                    y_resid_polylog = n_p(y_coef,abs(v_y()))-y_mean_flux();
                    subplot(2,2,4)
                    hold on
    %                 plot(v_y(y_mask),y_resid_polylog./y_mean_flux(y_mask),'ro')
    %                 plot(v_y(y_mask),(gfun(y_coef_gauss,(v_y(y_mask)))-y_mean_flux(y_mask))./y_mean_flux(y_mask),'bx')
                    plot(v_y,y_poly_rel_err,'ro')
                    plot(v_y,(gfun(y_coef_gauss,(v_y))-y_mean_flux)./y_mean_flux,'bx')
                    ylim([-1,1])
                    xlabel('$v_y$')
                    ylabel('$\Delta$')
                    title('$v_y$ residuals')

                    suptitle(sprintf('Single shot fit. T$_P$ =(%.2e,%.2e), T$_G$ =(%.2e,%.2e)',x_coef(3),y_coef(3),x_temp_gauss,y_temp_gauss));
                end % single shot plot
                shot_data.ok_mask(i) = 1;
            end %if enough counts



        catch
            warning('on')
            warning('Error in shot %u',i);
            shot_data.ok_mask(i) = 0;
        end %try shot
    %     cli_header(3,'Done.');

    %     end
    end %loop over shots

    %%
    % for these temperatures, k_B*T ~ 5e-30, and hbar*omega_trap ~ 1.5e-31
    % so this is right about the crossover
    ok_mask = logical(shot_data.ok_mask);
    idxs = logical(shot_data.ok_mask);
    colours = plasma(8);

    X_temps = 1e6*shot_data.temp(logical(ok_mask),1);
    Y_temps = 1e6*shot_data.temp(logical(ok_mask),2);
    X_g_temps = 1e6*shot_data.gauss_temp(logical(ok_mask),1);
    Y_g_temps = 1e6*shot_data.gauss_temp(logical(ok_mask),2);

    shotnums = shot_data.shot_num(logical(shot_data.ok_mask));
    shot_counts = shot_data.num_counts(ok_mask);
    % Pull temperature out of fit coefficient?
    % 1/(m * w_ho * poly_fac^(1/3)) = lambda_T 
    % w_ho = 2*pi*geomean([425,425,45]);
    % poly_fac = shot_data.amp/.08;
    % poly_fac = poly_fac(poly_fac < Inf & poly_fac > 0);
    % nanmean(poly_fac)
    % denom = nanmean((1./poly_fac).^(1/3));
    % lambda_T = denom/(mhe*w_ho)

    % % 


    cli_header(1,'Plotting...');
    stfig('PAL temp fits');
    clf
    subplot(2,2,1)
    hold on
    px=plot(shot_data.shot_num(logical(ok_mask)),X_temps,'.','MarkerSize',7,'Color',colours(1,:));
    py=plot(shot_data.shot_num(logical(ok_mask)),Y_temps,'.','MarkerSize',7,'Color',colours(4,:));
    gx=plot(shot_data.shot_num(logical(ok_mask)),X_g_temps,'.','MarkerSize',7,'Color',colours(5,:));
    gy=plot(shot_data.shot_num(logical(ok_mask)),Y_g_temps,'.','MarkerSize',7,'Color',colours(7,:));
    ylim([0,max(max([X_temps,Y_temps]))])
    xlabel('Shot number')
    ylabel('Temperature ($\mu$K)')
    % set(gca,'Yscale','log')
    ylim([0,3])
    legend([px,py,gx,gy],'polylog X','polylog Y','Gaussian X','Gaussian Y','Location','NorthEast')
    title('Fitted temps ($\mu$K)')

    % subplot(2,2,2)
    % hold on
    % histogram((shot_data.temp(:,1)-shot_data.temp(:,2))./nanmean(shot_data.temp,2),30,'FaceColor',[.1,.1,.7],'FaceAlpha',0.4)
    % xlabel('Relative err (x-y) / 2(x+y)')
    % ylabel('Num shots')
    % title('Relative difference in X,Y temps')


    subplot(4,2,2)
    hold on
    temp_scale = 1e6;
    px=plot(shot_data.num_counts(ok_mask),temp_scale*shot_data.temp(ok_mask,1),'.','MarkerSize',7,'Color',colours(1,:));
    gx=plot(shot_data.num_counts(ok_mask),temp_scale*shot_data.gauss_temp(ok_mask,1),'.','MarkerSize',7,'Color',colours(6,:));
    % legend([px,py,gx,gy],'polylog X','polylog Y','Gaussian X','Gaussian Y')
    legend('Polylog','Gaussian')
    ylim([0,1])
    % xlabel('Num counts')
    % xticks([])
    ylabel('X Temperature ($\mu$K)')
    title('Temperature-Number relations')

    subplot(4,2,4)
    hold on
    temp_scale = 1e6;
    py=plot(shot_data.num_counts(ok_mask),temp_scale*shot_data.temp(ok_mask,2),'.','MarkerSize',7,'Color',colours(1,:));
    gy=plot(shot_data.num_counts(ok_mask),temp_scale*shot_data.gauss_temp(ok_mask,2),'.','MarkerSize',7,'Color',colours(6,:));
    % legend([px,py,gx,gy],'polylog X','polylog Y','Gaussian X','Gaussian Y')
    % legend('Polylog','Gaussian')
    % set(gca,'Yscale','log')
    ylim([0.1,3])
    xlabel('Num counts')
    ylabel('Y Temperature ($\mu$K)')


    subplot(2,2,3)
    hold on
    px=plot(shotnums,1-shot_data.polylog_r_squared(ok_mask,1),'.','MarkerSize',7,'Color',colours(1,:));
    py=plot(shotnums,1-shot_data.polylog_r_squared(ok_mask,2),'.','MarkerSize',7,'Color',colours(4,:));
    gx=plot(shotnums,1-shot_data.gauss_r_squared(ok_mask,1),'.','MarkerSize',7,'Color',colours(6,:));
    gy=plot(shotnums,1-shot_data.gauss_r_squared(ok_mask,2),'.','MarkerSize',7,'Color',colours(7,:));
    set(gca,'Yscale','log')
    xlabel('Shot number')
    ylabel('$r^2$')
    title('Fit MSE')
    legend([px,py,gx,gy],'polylog X','polylog Y','Gaussian X','Gaussian Y')


    subplot(4,2,6)
    hold on
    temp_scale = 1e6;
    py=plot(shot_counts,shot_data.temp(ok_mask,1)-shot_data.temp(ok_mask,2),'.','MarkerSize',7,'Color',colours(1,:));
    py=plot(shot_counts,shot_data.gauss_temp(ok_mask,1)-shot_data.gauss_temp(ok_mask,2),'.','MarkerSize',7,'Color',colours(5,:));
    legend('Polylog','Gaussian')
    title('Agreement between XY fits')
    ylim([-1e-6,3e-6])
    xlabel('Shot number')
    ylabel('X-Y temp')

    subplot(4,2,8)
    hold on
    temp_scale = 1e6;
    py=plot(shot_counts,shot_data.temp(ok_mask,1)-shot_data.gauss_temp(ok_mask,1),'.','MarkerSize',7,'Color',colours(1,:));
    py=plot(shot_counts,shot_data.temp(ok_mask,2)-shot_data.gauss_temp(ok_mask,2),'.','MarkerSize',7,'Color',colours(5,:));
    legend('X','Y')
    title('Agreement between fit types')
    ylim([-1e-6,3e-6]   )
    xlabel('Shot number')
    ylabel('poly-gauss temp')


    suptitle('Pulse-aligned thermometry')

end