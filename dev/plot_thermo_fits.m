function plot_thermo_fits(data,opts)

    x_ok_mask = logical(data.temp.ok_mask(:,1));
    y_ok_mask = logical(data.temp.ok_mask(:,2));
    all_ok_mask = x_ok_mask & y_ok_mask;
    idxs = logical(data.temp.ok_mask);
    colours = plasma(8);
    temp_scale = 1e6;

    shotnums = data.temp.shot_num;
    shot_counts = data.temp.num_counts;

    % max_bose = 
    % ymax = 

    f=stfig('PAL temp fits');
    clf
    subplot(4,1,1)
    hold on
    px=plot(shotnums(all_ok_mask),temp_scale*data.temp.bose.temp(all_ok_mask,1),'.','MarkerSize',7,'Color',colours(1,:));
    py=plot(shotnums(all_ok_mask),temp_scale*data.temp.bose.temp(all_ok_mask,2),'.','MarkerSize',7,'Color',colours(4,:));
    gx=plot(shotnums(all_ok_mask),temp_scale*data.temp.gauss.temp(all_ok_mask,1),'.','MarkerSize',7,'Color',colours(5,:));
    gy=plot(shotnums(all_ok_mask),temp_scale*data.temp.gauss.temp(all_ok_mask,2),'.','MarkerSize',7,'Color',colours(7,:));
    %     ylim([0,max(max([X_temps,Y_temps]))])
    xlabel('Shot number')
    ylabel('Temperature ($\mu$K)')
    % set(gca,'Yscale','log')
    ylim([0,1])
    legend([px,py,gx,gy],'Bose X','Bose Y','Gaussian X','Gaussian Y','Location','NorthEast')
    title('Fitted temps ($\mu$K)')

    subplot(4,1,2)
    hold on
%     px=plot(shotnums(all_ok_mask),1-data.temp.bose.r_squared(all_ok_mask,1),'.','MarkerSize',7,'Color',colours(1,:));
%     py=plot(shotnums(all_ok_mask),1-data.temp.bose.r_squared(all_ok_mask,2),'.','MarkerSize',7,'Color',colours(4,:));
%     gx=plot(shotnums(all_ok_mask),1-data.temp.gauss.r_squared(all_ok_mask,1),'.','MarkerSize',7,'Color',colours(6,:));
%     gy=plot(shotnums(all_ok_mask),1-data.temp.gauss.r_squared(all_ok_mask,2),'.','MarkerSize',7,'Color',colours(7,:));
%     set(gca,'Yscale','log')
    plot(shotnums,shot_counts,'.','Color',colours(2,:))
    plot(shotnums(~all_ok_mask),shot_counts(~all_ok_mask),'r*')
    legend('Raw counts','Excluded shots')
    xlabel('Shot number')
%     ylabel('$1-r^2$')
    ylabel('Number of counts')
    title('MCP count trace')
%     legend([px,py,gx,gy],'Bose X','Bose Y','Gaussian X','Gaussian Y')

    % subplot(2,2,2)
    % hold on
    % histogram((data.temp.temp(:,1)-data.temp.temp(:,2))./nanmean(data.temp.temp,2),30,'FaceColor',[.1,.1,.7],'FaceAlpha',0.4)
    % xlabel('Relative err (x-y) / 2(x+y)')
    % ylabel('Num shots')
    % title('Relative difference in X,Y temps')

    subplot(2,2,3)
    hold on
    px=plot(data.temp.num_counts(all_ok_mask),temp_scale*data.temp.bose.temp(all_ok_mask,1),'.','MarkerSize',7,'Color',colours(1,:));
    gx=plot(data.temp.num_counts(all_ok_mask),temp_scale*data.temp.gauss.temp(all_ok_mask,1),'.','MarkerSize',7,'Color',colours(4,:));
    % legend([px,py,gx,gy],'Bose X','Bose Y','Gaussian X','Gaussian Y')
    legend([px,gx],'Bose X','Gaussian X')
    xlim([0,7e4])
    ylim([0,1])
    ylabel('Temperature ($\mu$K)')
    title('Temperature-Number correlation in X')

    subplot(2,2,4)
    hold on
    py=plot(shot_counts(all_ok_mask),temp_scale*data.temp.bose.temp(all_ok_mask,2),'.','MarkerSize',7,'Color',colours(6,:));
    gy=plot(shot_counts(all_ok_mask),temp_scale*data.temp.gauss.temp(all_ok_mask,2),'.','MarkerSize',7,'Color',colours(7,:));
    legend([py,gy],'Bose Y','Gaussian Y')
    % legend('Polylog','Gaussian')
    % set(gca,'Yscale','log')
    xlim([0,7e4])
    ylim([0,1])
    xlabel('Num counts')
    ylabel('Temperature ($\mu$K)')
    title('Temperature-Number correlation in Y')
    
    suptitle(sprintf('Temperature fits: $%s$',opts.plt_label))
    
    filename = sprintf('temp_fits_%s',opts.plt_label);
    saveas(f,[filename,'.fig']);
    saveas(f,[filename,'.png']);
    saveas(f,[filename,'.svg']);
    saveas(f,[filename,'.epsc']);
end


% 
% 
% subplot(2,2,4)
% hold on
% plot(shot_counts(all_ok_mask),data.temp.bose.temp(all_ok_mask,1)-data.temp.bose.temp(all_ok_mask,2),'.','MarkerSize',7,'Color',colours(1,:));
% plot(shot_counts(all_ok_mask),data.temp.gauss.temp(all_ok_mask,1)-data.temp.gauss.temp(all_ok_mask,2),'.','MarkerSize',7,'Color',colours(4,:));
% legend('Polylog','Gaussian')
% title('Agreement between XY fits')
% ylim([-1e-6,3e-6])
% xlabel('Shot number')
% ylabel('X-Y temp')
% 
% plot(shot_counts(x_ok_mask),data.temp.bose.temp(x_ok_mask,1)-data.temp.gauss.temp(x_ok_mask,1),'.','MarkerSize',7,'Color',colours(6,:));
% plot(shot_counts(y_ok_mask),data.temp.bose.temp(y_ok_mask,2)-data.temp.gauss.temp(y_ok_mask,2),'.','MarkerSize',7,'Color',colours(7,:));
% legend('X','Y')
% title('Agreement between fit types')
% ylim([-1e-6,3e-6]   )
% xlabel('Shot number')
% ylabel('Bose-gauss temp')
% suptitle('Pulse-aligned thermometry')