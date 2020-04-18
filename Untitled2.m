cli_header(2,'Plotting:');
        stfig('Temperature fits');
        clf;
        subplot(4,2,1)
        hold on
        plot(keep_idx(x_mask), temps(x_mask,1),'k:')
        plot(keep_idx(y_mask), temps(y_mask,2),'r:')
        box off
        xlabel('Shot number')
        ylim([0,1.5*max(max(temps(net_mask)))])
        xticks([])
        ylabel('Fit width')
        title('Time series')

        subplot(4,2,3)
        n_plot = data_tdc.num_counts(x_mask | y_mask);
        plot(data_tdc.shot_num(x_mask | y_mask),n_plot-max(n_plot),'b:')
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
        cli_header(2,'Done.');