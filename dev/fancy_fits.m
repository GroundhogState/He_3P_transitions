function data = fancy_fits(data,opts)
cli_header({0,'Fitting peaks...'})
num_cats = numel(data.cat);
for cidx = 1:num_cats
    opts.fig_idx = cidx;
    data.cat{cidx}.pfits = fancy_fits_core(data.cat{cidx},opts);
end

cli_header({1,'Done.'})
end


function pfits = fancy_fits_core(data,opts)

    if isfield(opts,'fig_idx')
        fnum = 54714 + opts.fig_idx;
    else
        fnum = 54714;
    end
    signal = data.spec.signal;
    sfigure(fnum+10);
    plot(signal);
    hold on
    
    
    f1=sfigure(fnum);
    clf;
    subplot(2,1,1)
    plot(data.spec.freq,data.spec.signal,'kx')
    hold on
    xlabel('Frequency (MHz)')
    ylabel('Number loss')
    title(['Fitting peaks in ',strrep(data.data.class{1},'_',' ')])
    
    npeaks = length(data.peaks.peak_num);
%     data.cat{cidx}.pfits = zeros(npeaks,1);
    [~,peak_order] = sort(data.peaks.prominences,'descend');
    

    
    for idx = 1:npeaks
        pidx = peak_order(idx);
        
        f_cen = data.peaks.freqs(pidx);
    
        f_shifted = data.spec.freq - f_cen;
        f_win_mask = f_shifted < opts.peak.fitwidth & f_shifted > -opts.peak.fitwidth;
        f_fit = f_shifted(f_win_mask);
        sig_fit = signal(f_win_mask);
        
        % Fit the data
        f_disp = linspace(min(f_fit),max(f_fit),100);
        f_guesses = data.peaks.freqs(pidx)-f_cen;
        h_guesses = data.peaks.vals(pidx);
        w_guesses = data.peaks.widths(pidx)/2;
        gfun = @(p,x) p(1)*exp(-((x-p(2))/p(3)).^2);
        lfun = @(p,x) p(1)*p(3)./((x-p(2)).^2 + (p(3)).^2);
        g0 = [h_guesses,f_guesses,sqrt(w_guesses)];
        l0 = [0.5*w_guesses*(1-min(data.spec.signal)),f_guesses,0.5*w_guesses];
        gfit = fitnlm(f_fit ,sig_fit,gfun,g0);
        lfit = fitnlm(f_fit ,sig_fit,lfun,l0);
        lf_plot=lfun(lfit.Coefficients.Estimate,f_disp );

        
        % Write the output
        pfits.lorz{pidx} = lfit;
        pfits.gaus{pidx} = gfit;
        pfits.offset(pidx) = f_cen;
        pfits.lor_prms(pidx,:) = [lfit.Coefficients.Estimate(2)+f_cen,2*lfit.Coefficients.Estimate(3)]; %cen, width
        pfits.lor_err(pidx,:) = [lfit.Coefficients.SE(2),2*lfit.Coefficients.SE(3)]; %cen, width

        
        % Plot the results
        sfigure(fnum);
        subplot(2,1,1)
        plot(f_disp+f_cen ,gfun(gfit.Coefficients.Estimate,f_disp ),'b-')
        plot(f_disp+f_cen ,lf_plot,'r-')
        legend('raw data','Gaussian fit','Lorentzian fit')
        


        subplot(2,1,2)
        plot(f_fit+f_cen,gfun(gfit.Coefficients.Estimate,f_fit)-sig_fit,'bx')
        hold on
        plot(f_fit+f_cen,lfun(lfit.Coefficients.Estimate,f_fit)-sig_fit,'rx')
        legend('Gaussian residuals','Lorentzian residuals')
        xlabel('Frequency (MHz)')
        title('Fit error')
        
        % Print output
        fprintf(' --- Peak from stage %u\n',opts.fig_idx)
        fprintf('Peak centre        %.2f(%.2f) MHz\n',lfit.Coefficients.Estimate(2)+f_cen,lfit.Coefficients.SE(2));
        fprintf('Peak width         %.2f(%.2f) MHZ\n',2*lfit.Coefficients.Estimate(3),2*lfit.Coefficients.SE(3));%p(3)=0.5*Gamma
        fprintf('Peak height        %.2f(%.2f)\n',lfit.Coefficients.Estimate(1),2*lfit.Coefficients.SE(1));
    
        sfigure(fnum+10);
        signal(f_win_mask) = signal(f_win_mask) - lfun(lfit.Coefficients.Estimate,f_shifted(f_win_mask));
        plot(signal);

    end
    % Save the figures
    imname = sprintf('peak_fit_%u',opts.fig_idx);
    filename1 = fullfile(opts.out_dir,imname);
    saveas(f1,[filename1,'.fig']);
    saveas(f1,[filename1,'.png'])
    
end
