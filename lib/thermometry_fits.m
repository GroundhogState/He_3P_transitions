function shot_data = thermometry_fits(data,opts)

    cache_opts=[];
    if isfield(opts.temp,'cache_opts'), cache_opts=opts.temp.cache_opts; end
    cache_opts.verbose=1;

    opts.shot_range = nan;
    
    
    opts.shotnums = data.tdc.shot_num;
    if ~isnan(opts.temp.shot_nums)
        shot_nums = opts.temp.shot_nums;
        shot_nums = shot_nums(shot_nums < max(opts.shotnums) & shot_nums > 1);
        opts.shotnums = shot_nums;        
    end
    cli_header(1,'Fitting thermal profiles');
    shot_data = thermo_fit_core(data,opts);
    

    %%
    % for these temperatures, k_B*T ~ 5e-30, and hbar*omega_trap ~ 1.5e-31
    % so this is right about the crossover

    cli_header(2,'Done.');
end

function shot_data = thermo_fit_core(data,opts)
    mcp_data = data.tdc;
    mcp_data.all_ok = ones(size(mcp_data.N_atoms));
    num_shots = length(opts.shotnums);
    opts.global.fall_velocity = 9.796*opts.const.fall_time;
    opts.global.fall_time = opts.const.fall_time;

    % Options for binning and show_TXY
    opts.t0 = 0.4145; % Centre of first pulse (sec)
    opts.pulsedt = .008; %time between pulses (sec)
    opts.num_bins = [50,50,50];
    opts.hist_2d = false;
    opts.draw_plots = false;
    opts.centre_bec = true;
    opts.verbose = false;
    opts.xylim = [-0.03,.03;-.03,.03];
    opts.pulse_twindow = .003; 
    
    min_flux = 1e4;
    num_pulses = 150;
    single_shot_plot = opts.single_shot_plot;
    do_refits = true;
    min_shot_counts = 1e2;

    x_refit_lvl = 0.15;
    y_refit_lvl = 0.15;

    const.g0 = 9.796;
    const.mhe = 6.6433e-27;
    const.kb=1.380e-23;
    
    A = 1e3;
    T_guess=2e-7;
    fit_guess = [A,0,T_guess,0];
    g_tol = 1e-5;
    n_p = @(p,v) p(1)*g_bose(exp(-(const.mhe*(v-p(2))).^2/(2*const.mhe*const.kb*p(3))),g_tol)+p(4);
    gfun = @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3)))+b(4); 
    fit_options = statset('TolFun',1e-6);
    
    shot_data=[];
    shot_data.pulse_data = cell(num_shots,1);
    shot_data.ok_mask = zeros(num_shots,2);
    shot_data.shot_num= ones(num_shots,1);
    shot_data.num_counts= ones(num_shots,1);
    
    for idx=1:num_shots
        i = opts.shotnums(idx);
        warning('off')
        if mod(idx-1,10)==0
            cli_header('Working on shot %u/%u...',idx,num_shots);
        end
        try
            mask = zeros(size(data.tdc.N_atoms));
            mask(i) = 1;
            tshot = struct_mask(data.tdc,logical(mask));
            if tshot.N_atoms < min_shot_counts
                warning('on')
                warning('Insufficient counts in shot %u',i)
                warning('off')
                shot_data.ok_mask(i,:) = 0;
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

                X = pulse_data.x_cens(1,:);
                x_mean_flux = mean(pulse_data.x_flux);
                v_x = X/opts.const.fall_time;

                Y = pulse_data.y_cens(1,:);
                y_mean_flux = mean(pulse_data.y_flux);
                v_y = Y/opts.const.fall_time; 

                x_mask = abs(v_x) > .01;
                if sum(x_mask) < 10
                    warning('on')
                    warning('Too few points to fix in shot %u X',i)
                    shot_data.ok_mask(i,1) = 0;
                else
                    x_mdl = fitnlm(v_x(x_mask),x_mean_flux(x_mask),n_p,fit_guess,...
                                'Options',fit_options);
                    x_coef.bose= x_mdl.Coefficients.Estimate;
                    x_SE = x_mdl.Coefficients.SE;
                    bose_resid_x =(n_p(x_coef.bose,v_x)-x_mean_flux);
                    x.bose.rel_err = bose_resid_x./x_mean_flux;

                     % % Fit Gaussian functions
                    x_guess.gauss=[mean(x_mean_flux),0,.0007,0];
                    x_mdlfit.gauss = fitnlm(v_x(x_mask),x_mean_flux(x_mask),gfun,x_guess.gauss,...
                                'Options',fit_options);
                    x_coef.gauss = x_mdlfit.gauss.Coefficients.Estimate;
                    gauss_resid_x = (gfun(x_coef.gauss,v_x)-x_mean_flux);
                    x.gauss.rel_err = gauss_resid_x./x_mean_flux;   

                    %re-fit
                    x_refit_mask = abs(x.gauss.rel_err) < x_refit_lvl | abs(x.bose.rel_err) < x_refit_lvl;
                    if sum(x_refit_mask) > sum(x_mask) && do_refits
                        % New things to fit
                        x_mask = x_refit_mask & abs(v_x) > .01;
                        x_mdl = fitnlm(v_x(x_mask),x_mean_flux(x_mask),n_p,fit_guess,...
                                    'Options',fit_options);
                        x_coef.bose= x_mdl.Coefficients.Estimate;
                        x_SE = x_mdl.Coefficients.SE;
                        bose_resid_x =(n_p(x_coef.bose,v_x)-x_mean_flux);
                        x.bose.rel_err = bose_resid_x./x_mean_flux;

                         % % Fit Gaussian functions
                        x_guess.gauss=[mean(x_mean_flux),0,.0007,0];
                        x_mdlfit.gauss = fitnlm(v_x(x_mask),x_mean_flux(x_mask),gfun,x_guess.gauss,...
                                    'Options',fit_options);
                        x_coef.gauss = x_mdlfit.gauss.Coefficients.Estimate;
                        gauss_resid_x = (gfun(x_coef.gauss,v_x)-x_mean_flux);
                        x.gauss.rel_err = gauss_resid_x./x_mean_flux; 
                    end
                    shot_data.ok_mask(i,1) = 1;
                end
                
                % Fitting Y
                % Polylog
                y_mask = abs(v_y) > 0.034 & y_mean_flux > 50;
                if sum(y_mask) < 10
                    warning('on')
                    warning(sprintf('Too few points to fix in shot %u Y',i))
                    shot_data.ok_mask(i,2) = 0;
                    warning('off')
                else
                    y_mdl = fitnlm(v_y(y_mask),y_mean_flux(y_mask),n_p,fit_guess,...
                            'Options',fit_options);
                    y_coef.bose = y_mdl.Coefficients.Estimate;
                    y_SE = y_mdl.Coefficients.SE;
                    bose_resid_y = (n_p(y_coef.bose,v_y)-y_mean_flux);
                    y.bose.rel_err = bose_resid_y./y_mean_flux;       
                    % Gaussian
                    y_guess_guess=[mean(x_mean_flux),0,.0007,0];
                    y_mdlfit.gauss = fitnlm(v_y(y_mask),y_mean_flux(y_mask),gfun,y_guess_guess,...
                            'Options',fit_options);
                    y_coef.gauss = y_mdlfit.gauss.Coefficients.Estimate;
                    gauss_resid_y = (gfun(y_coef.gauss,v_y)-y_mean_flux);
                    y.gauss.rel_err = gauss_resid_y./y_mean_flux; 

                    %re-fit
                    y_refit_mask = (abs(y.gauss.rel_err) < y_refit_lvl | abs(y.bose.rel_err) < y_refit_lvl) & abs(v_y) > 0.01;
                    if sum(y_refit_mask) > sum(y_mask) && do_refits
                        % New things to fit
                        y_mask = y_refit_mask;
                        y_mdl = fitnlm(v_y(y_mask),y_mean_flux(y_mask),n_p,fit_guess,...
                                    'Options',fit_options);
                        y_coef.bose = y_mdl.Coefficients.Estimate;
                        y_SE = y_mdl.Coefficients.SE;
                        bose_resid_y = (n_p(y_coef.bose,v_y)-y_mean_flux);
                        y.bose.rel_err = bose_resid_y./y_mean_flux;
                        % Gaussian
                        y_guess_guess=[mean(x_mean_flux),0,.0007,0];
                        y_mdlfit.gauss = fitnlm(v_y(y_mask),y_mean_flux(y_mask),gfun,y_guess_guess,...
                                    'Options',fit_options);
                        y_coef.gauss = y_mdlfit.gauss.Coefficients.Estimate;
                        gauss_resid_y = (gfun(y_coef.gauss,v_y)-y_mean_flux);
                        y.gauss.rel_err = gauss_resid_y./y_mean_flux; 
                    end
                    shot_data.ok_mask(i,2) = 1;
                end

                bose.r_squared = [x_mdl.Rsquared.Ordinary,y_mdl.Rsquared.Ordinary];
                gauss.r_squared = [x_mdlfit.gauss.Rsquared.Ordinary,y_mdlfit.gauss.Rsquared.Ordinary];
                gauss_temp_x =(abs(x_coef.gauss(3))) *const.mhe/const.kb;
                gauss_temp_y = (abs(y_coef.gauss(3))) *const.mhe/const.kb;

    %             [~, msgid] = lastwarn; %were there warnings in that attempt?
    %             if strcmp(msgid,'stats:nlinfit:IllConditionedJacobian')
    %                 cli_header(2,'Shot %u Gaussian poorly fitted\n',i);
    %     %             shot_data.ok_mask(i) = 0;
    %                 warning('Restting warning...');
    %             end


                pulse_data.gauss.residual.x = gauss_resid_x;
                pulse_data.gauss.residual.y = gauss_resid_y;
                pulse_data.bose.residual.x = bose_resid_x;
                pulse_data.bose.residual.y = bose_resid_y;
                pulse_data.x_mask = x_mask;
                pulse_data.y_mask = y_mask;
                
                shot_data.pulse_data{i} = pulse_data;
                shot_data.shot_num(i) = i;
                shot_data.num_counts(i) = mcp_data.N_atoms(i);
                
                shot_data.bose.temp(i,:) = [x_coef.bose(3),y_coef.bose(3)];
                shot_data.bose.amp(i,:) = [x_coef.bose(1),y_coef.bose(1)];
                shot_data.bose.r_squared(i,:) = bose.r_squared;
                
                shot_data.gauss.temp(i,:) = [gauss_temp_x,gauss_temp_y];
                shot_data.gauss.amp(i,:) = [x_coef.gauss(1),y_coef.gauss(1)];
                shot_data.gauss.r_squared(i,:) = gauss.r_squared;
                
                if single_shot_plot   
                    cli_header(2,'Plotting...');
                    stfig('Pulse data');
                    clf
                    subplot(2,2,1)
                    hold on
                    plot(v_x,x_mean_flux,'k:')
                    plot(v_x(x_mask),x_mean_flux(x_mask),'k*')
                    plot(v_x(x_refit_mask),x_mean_flux(x_refit_mask),'ko')
                    plot(v_x,n_p(x_coef.bose,abs(v_x)),'r')
                    plot(v_x,gfun(x_coef.gauss,v_x),'b')
                    ylim([100,max(2*x_mean_flux)])
                    legend('Flux','Fit domain 1','Fit domain 2','Polylog','Gaussian')
                    xlabel('$v_x$')
                    ylabel('Flux')
                    set(gca,'Yscale','log')
                    title('Mean X profile')

                    subplot(2,2,2)
                    hold on
                    plot(v_y,y_mean_flux,'k:')
                    plot(v_y(y_mask),y_mean_flux(y_mask),'k*')
                    plot(v_y(y_refit_mask),y_mean_flux(y_refit_mask),'ko')
                    plot(v_y,n_p(y_coef.bose,abs(v_y)),'r')
                    plot(v_y,gfun(y_coef.gauss,v_y),'b')
                    xlabel('$v_y$')
                    ylabel('Flux')
                    set(gca,'Yscale','log')
                    title('Mean Y profile')

                    subplot(2,2,3)
                    hold on
                    plot(v_x,(n_p(x_coef.bose,abs(v_x))-x_mean_flux)./x_mean_flux,'ro')
                    plot(v_x,(gfun(x_coef.gauss,(v_x))-x_mean_flux)./x_mean_flux,'bx')
                    legend('Polylog','Gaussian','Location','NorthEast')
                    ylim([-1,1])
                    xlabel('$v_x$')
                    ylabel('$\Delta$')
                    title('$v_x$ residuals')

                    subplot(2,2,4)
                    hold on
                    plot(v_y,y.bose.rel_err,'ro')
                    plot(v_y,(gfun(y_coef.gauss,(v_y))-y_mean_flux)./y_mean_flux,'bx')
                    ylim([-1,1])
                    xlabel('$v_y$')
                    ylabel('$\Delta$')
                    title('$v_y$ residuals')

                    suptitle(sprintf('Shot %u: T$_P$ =(%.2e,%.2e), T$_G$ =(%.2e,%.2e), N=%u',i,x_coef.bose(3),y_coef.bose(3),gauss_temp_x,gauss_temp_y,tshot.N_atoms));
                end % single shot plot
                
            end %if enough counts

        catch
            warning('on')
            warning('Error in shot %u',i);
            shot_data.ok_mask(i,:) = 0;
        end %try shot
    end %loop over shots
end