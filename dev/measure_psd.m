function psd_data = measure_psd(data_tdc,opts)
% A function that fits the temperature profiles of PAL pulses. 
% Needs: t_start, t_binwidth, num_peaks
%   or?  t_start, t_high, t_low, num_peaks
  % or!  Start time & profile of RF -

cli_header(1,'Fitting temperatures...');
    lims = [0.38,2.1;
            -.03,.03;
            -.03,.03];
    min_num_counts = 1e4;
    count_mask = data_tdc.N_atoms > min_num_counts;
    keep_data = data_tdc.counts_txy;
    keep_counts = data_tdc.num_counts;
    
    num_bins = [200,200,200];
    threshold = 90;
    
    keep_idx = data_tdc.shot_num;
    num_shots = length(keep_data);
    if ~isnan(opts.shots_chosen)
        shots_chosen = opts.shots_chosen;
    else
        shots_chosen = 1:num_shots;
    end

    num_shots_chosen = length(shots_chosen);

    gfun = @(p,x) p(1)*exp(-((x-p(2))/p(3)).^2)+p(4);
    temps = nan(num_shots,2);
    temps_SE = nan(num_shots,2);

    not_ok_mask = zeros(num_shots_chosen,1);
    
    stfig('Temperature fits');
    clf;
%     profile on
    for idx = 1:num_shots_chosen
        test_idx = shots_chosen(idx);
        if mod(idx,10) == 0
            cli_header(2,'Shot # %u/%u',idx,num_shots_chosen);
%         fprintf('Shot number %u:\n',shots_chosen(idx))
        end
        if count_mask(test_idx)
            try
            these_counts = keep_data{test_idx};
            these_counts = masktxy_square(these_counts,lims);
            these_counts = these_counts - mean(these_counts);

            x_data = these_counts(:,2);
            [x_counts,x_edges] = histcounts(x_data,num_bins(2));
            x_cens = 0.5*(x_edges(2:end) + x_edges(1:end-1));
            x_mask = x_counts < threshold;
            x_guess = [1.5*mean(x_counts),0,2*std(x_data),0];
            x_fit = fitnlm(x_cens(x_mask),x_counts(x_mask),gfun,x_guess);

            y_data = these_counts(:,3);
            [y_counts,y_edges] = histcounts(y_data,num_bins(3));
            y_cens = 0.5*(y_edges(2:end) + y_edges(1:end-1));
            y_mask = y_counts < threshold;
            y_guess = [1.5*mean(y_counts),0,2*std(y_data),0];
            y_fit = fitnlm(y_cens(y_mask),y_counts(y_mask),gfun,y_guess);

            temps(test_idx,:) = [x_fit.Coefficients.Estimate(3),y_fit.Coefficients.Estimate(3)];
            temps_SE(test_idx,:) = [x_fit.Coefficients.SE(3),y_fit.Coefficients.SE(3)];
            catch
                fprintf('Fit failed in shot %u\n',test_idx)        
                not_ok_mask(idx) = 1;
            end %try fits
            if opts.single_plot && idx == 1
                x_binsize = diff(x_edges);
                y_binsize = diff(y_edges);

                cli_header(2,'Plotting shot %u ...',test_idx)
                stfig('PAL Profiles');
                clf

                subplot(2,2,1)
                hold on
                plot(x_cens(x_mask),x_counts(x_mask)./x_binsize(x_mask),'b')
                plot(y_cens(y_mask),y_counts(y_mask)./y_binsize(y_mask),'k')

                subplot(2,2,2)
                hold on
                plot(x_cens(x_mask),x_counts(x_mask),'b')
                plot(y_cens(y_mask),y_counts(y_mask),'k')

                subplot(2,2,4)
                plot(x_cens,x_counts./y_counts)

                subplot(2,2,3)
                hold on
                plot(x_cens,x_counts)
                plot(y_cens,y_counts)
                plot(x_cens,gfun(x_guess,x_cens),'k:')
                plot(y_cens,gfun(y_guess,y_cens),'r:')
                plot(x_cens,gfun(x_fit.Coefficients.Estimate,x_cens),'k')
                plot(y_cens,gfun(y_fit.Coefficients.Estimate,y_cens),'r')
                set(gca,'Yscale','log')

                subplot(2,2,4)
                hold on
                suptitle(sprintf('Example fit, %u counts',keep_counts(test_idx)))
            end       % single_plot graphics
                
            
        else %if num too low
            cli_header(2,'Skipping shot %u, num too low',test_idx);
        end  % num check
    end %loop over shots

%     count_mask = data_tdc.num_counts > 2e4;

%     master_mask = ~not_ok_mask;
    x_physical_mask = temps_SE(:,1) < 1e-3  & temps(:,1) > 0 & temps(:,1) < 1;
    x_relative_unc = temps_SE(:,1)./temps(:,1);
    x_unc_mask = x_relative_unc<.05;
    x_mask = x_unc_mask & x_physical_mask;

    y_physical_mask = temps_SE(:,2) < 1e-3  & temps(:,2) > 0 & temps(:,2) < 1;
    y_relative_unc = temps_SE(:,2)./temps(:,2);
    y_unc_mask = y_relative_unc<.05;
    y_mask = y_unc_mask & y_physical_mask;

    net_mask = x_mask & y_mask & count_mask;
 
%     master_mask = ~not_ok_mask;


    if num_shots_chosen>5
        
        cli_header(2,'Plotting:');
        stfig('Temperature fits');
        clf;
        subplot(4,2,1)
        % yyaxis left
        hold on
        plot(keep_idx(x_mask), temps(x_mask,1),'k:')
        plot(keep_idx(y_mask), temps(y_mask,2),'r:')
        box off
        % xlabel('Shot number')
        ylim([0,2*max(max(temps))])
        xticks([])
        ylabel('Fit width')
        title('Time series')

        subplot(4,2,3)
        n_plot = data_tdc.num_counts(x_mask | y_mask);
        % yyaxis right
        plot(data_tdc.shot_num(x_mask | y_mask),n_plot-max(n_plot),'b.')
        box off
        xlabel('Shot number')
        ylabel('Atoms detected')

        subplot(2,2,2)
        hold on
        plot(keep_counts(x_mask & y_mask),x_relative_unc(x_mask & y_mask),'k.')
        plot(keep_counts(x_mask & y_mask),y_relative_unc(x_mask & y_mask),'r.')
        xlabel('Counts')
        ylabel('Relative error')
        legend('X fit','Y fit')
        title('Relative error')

        subplot(2,2,3)
        hold on
        plot(temps(net_mask,1),x_relative_unc(x_mask & y_mask),'k.')
        plot(temps(net_mask,2),y_relative_unc(x_mask & y_mask),'r.')
        xlabel('Fit width')
        ylabel('Absolute err uncertainty ')
        legend('X fit','Y fit')
        title('Absolute error')

        subplot(2,2,4)
        hold on
        plot(keep_counts(net_mask),temps(net_mask,1),'k.')
        plot(keep_counts(net_mask),temps(net_mask,2),'r.')
        xlabel('Atoms detected')
        ylabel('Fit width')
%         set(gca,'Yscale','log')
%         set(gca,'Xscale','log')
        title('Num-width  correlation')
        legend('X fit','Y fit')
        cli_header(2,'Done.')

        suptitle(sprintf('XY Gaussian fits: $5^3D_{2,3}$ lines'))
    end
%     keep_out = struct_mask(,count_mask);
%     psd_data = struct_mask(data_tdc);
    psd_data.temps = temps;
    psd_data.temps_SE= temps_SE;
    psd_data.net_mask = net_mask;
    
    cli_header(1,'Done.');
end



