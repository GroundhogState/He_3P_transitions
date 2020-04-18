function al_pulses=fit_gauss_to_al_pulses(anal_opts,mcp_data)
% warning al_pulses.pos_stat will be removed in a furture release
% global const

% Inputs (required):
%     mcp_data.counts_txy
%     mcp_data.all_ok
% 
%     anal_opts.pulsedt
%     anal_opts.pulses
%     anal_opts.t0
%     anal_opts.start_pulse
%     anal_opts.fit.hist_sigma
%     anal_opts.xylim
%     anal_opts.pulse_twindow
%     anal_opts.global.fall_velocity
%     anal_opts.global.fall_time
%     anal_opts.plot.all

iimax=size(mcp_data.counts_txy,2);
al_pulses=[];
al_pulses.pulsedt=anal_opts.pulsedt;
al_pulses.window=nan(anal_opts.pulses,3,2); %initalize
al_pulses.num_counts=nan(iimax,anal_opts.pulses);
al_pulses.pos.mean=nan(iimax,anal_opts.pulses,3);
al_pulses.pos.std=nan(iimax,anal_opts.pulses,3);

al_pulses.fit.Estimate=nan(iimax,anal_opts.pulses,4);
al_pulses.fit.SE=al_pulses.fit.Estimate;
al_pulses.fit.temperature.val=nan(iimax,anal_opts.pulses);
al_pulses.fit.temperature.unc=nan(iimax,anal_opts.pulses);

fig_single_pulse=stfig('single pulse fit');

%gauss_fun1d= amp*exp(-1*((x1-mu)^2)/(2*sig^2))
% gauss_fun1d = @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2));
% coeff_names={'amp','mu','sig'};
gauss_fun1d = @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2))+b(4); 
coeff_names={'amp','mu','sig','offset'};

fprintf('fitting A.L. pulses in files %04u:%04u',size(mcp_data.counts_txy,2),0)
first_good_shot=true;
for shot=2:iimax
        shot_txy_counts=mcp_data.counts_txy{shot};
        if mcp_data.all_ok(shot) & numel(shot_txy_counts)>1e3;
            for pulse=1:anal_opts.pulses
                %set up time window centered arround t0
                t_pulse_cen=anal_opts.t0+anal_opts.pulsedt...
                    *(anal_opts.start_pulse+pulse-2);
                trange=t_pulse_cen+anal_opts.pulse_twindow*[-0.5,0.5];
                pulse_win_txy=[trange;anal_opts.xylim]; 
                counts_pulse=masktxy_square(shot_txy_counts,pulse_win_txy);
                
                %now we make a smooth histogram in time and fit
                shist_out=smooth_hist(counts_pulse(:,1),'lims',trange,'sigma',anal_opts.fit.hist_sigma);
                
                xdata=shist_out.bin.centers-t_pulse_cen;
                ydata=shist_out.count_rate.smooth;
                if sum(ydata>0)>5
                    amp_guess=max(ydata);
                    mu_guess=wmean(xdata,ydata); %compute the weighted mean
                    sig_guess=sqrt(sum((xdata-mu_guess).^2.*ydata)/sum(ydata)); %compute the mean square weighted deviation
                    fo = statset('TolFun',10^-6,...
                        'TolX',1e-4,...
                        'MaxIter',1e4,...
                        'UseParallel',1);
                    % 'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',...
                    inital_guess=[amp_guess,mu_guess,sig_guess,0];
                    fitobject=fitnlm(xdata,ydata,...
                        gauss_fun1d,...
                         inital_guess,...
                        'CoefficientNames',coeff_names,'Options',fo);
                    fit_coeff=fitobject.Coefficients.Estimate;
                    fit_se=fitobject.Coefficients.SE;
                    al_pulses.fit.Estimate(shot,pulse,:)=fit_coeff;
                    al_pulses.fit.SE(shot,pulse,:)=fit_se;
                    al_pulses.fit.RMSE(shot,pulse)=fitobject.RMSE;
                    al_pulses.fit.x_offset=t_pulse_cen;

                    temperature_val=(abs(fit_coeff(3))*anal_opts.global.fall_velocity/anal_opts.global.fall_time)^2 *const.mhe/const.kb;
                    temperature_unc=temperature_val*2*fit_se(3)/abs(fit_coeff(3));

                    al_pulses.fit.temperature.val(shot,pulse)=temperature_val;
                    al_pulses.fit.temperature.unc(shot,pulse)=temperature_unc;

                    if anal_opts.plot.all
                        yscaling=1e-3;
                        stfig(fig_single_pulse);
                        clf
                        plot(shist_out.bin.centers,shist_out.count_rate.smooth*yscaling,'k')
                        xlabel('t, time(s)');
                        ylabel('Count Rate(kHz/File)');
                        x_sample_fit=col_vec(linspace(min(xdata),max(xdata),1e3));
                        [ysamp_val,ysamp_ci]=predict(fitobject,x_sample_fit);
                        x_sample_fit=x_sample_fit+t_pulse_cen;
                        hold on
                        plot(x_sample_fit,ysamp_ci*yscaling,'color',[1,1,1].*0.5)
                        plot(x_sample_fit,ysamp_val*yscaling,'r')
                        plot(x_sample_fit,gauss_fun1d(inital_guess,x_sample_fit-t_pulse_cen)*yscaling)


                        cen_units='ms';
                        cen_str=string_value_with_unc(fit_coeff(2)*1e3,fit_se(2)*1e3,'b');
                        width_units='ms';
                        width_str=string_value_with_unc(fit_coeff(3)*1e3,fit_se(3)*1e3,'b');
                        offset_units='kHz';
                        offset_str=string_value_with_unc(fit_coeff(4)*1e-3,fit_se(4)*1e-3,'b');

                        temperature_str=string_value_with_unc(1e6*temperature_val,1e6*temperature_unc,'b');
                        str=sprintf('Gauss fit \n   Cen %s %s \n   Width %s %s \n   offset %s %s\n   Temp.(no interactions)%s uk',...
                            cen_str,cen_units,width_str,width_units,offset_str,offset_units,temperature_str);
                        text(0.01,0.9,str,'Units','normalized'); 
                        hold off
                        drawnow

                    end
                end
                if first_good_shot
                    %only need to store this on first shot becasue the same for
                    %all shots
                    al_pulses.window(pulse,:,:)=pulse_win_txy; 
                    al_pulses.time_cen(pulse,:)=t_pulse_cen;
                end
                al_pulses.num_counts(shot,pulse)=size(counts_pulse(:,3),1);
                al_pulses.pos.mean(shot,pulse,:)=mean(counts_pulse,1);
                al_pulses.pos.std(shot,pulse,:)=std(counts_pulse,1);
                al_pulses.pos_stat(shot,pulse,:)=[...
                                           mean(counts_pulse(:,1)),...
                                           mean(counts_pulse(:,2)),...
                                           mean(counts_pulse(:,3)),...
                                           std(counts_pulse(:,1)),...
                                           std(counts_pulse(:,2)),...
                                           std(counts_pulse(:,3))]; 

            end%pulse
            if first_good_shot,first_good_shot=false; end
        %check that the calculated t0 is close enough
        
        tol_t0_match=1e-2; %in factors of bin size
        tol_t0_match=tol_t0_match*anal_opts.pulsedt;
        mean_cen_time=mean(al_pulses.pos.mean(shot,:,1)-al_pulses.time_cen');
        if abs(mean_cen_time)>tol_t0_match
            est_t0=anal_opts.t0+mean_cen_time;
            warning('pulses are not centered in time pehaps t0 should be %.5f',est_t0)
        end
        end%is mcp_data.all_ok
        if mod(shot,1)==0,fprintf('\b\b\b\b%04u',shot),end     
%to set the pulse t0 right it can be handy to uncomment the next line
%
end%shots

%check if the average difference between the mean count position (vs the expected center) over all shots is outside a tolerance
tol_t0_match=1e-3; %in factors of bin size
tol_t0_match=tol_t0_match*anal_opts.pulsedt; 
mean_cen_time=nanmean(arrayfun(@(shotnum) mean(al_pulses.pos.mean(shotnum,:,1)-al_pulses.time_cen'),1:size(al_pulses.pos.mean,1)));
if abs(mean_cen_time)>tol_t0_match
    est_t0=anal_opts.t0+mean_cen_time;
    fprintf('t0 should be %.6f',est_t0)
end

fprintf('...Done\n') 


end



