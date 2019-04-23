function data = fit_detected_peaks(data,opts)
header({0,'Fitting peaks...'})
num_cats = numel(data.cat);
for cidx = 1:num_cats
    opts.fig_idx = cidx;
    data.cat{cidx}.pfits = fit_found_peak_core(data.cat{cidx},opts);
end

header({1,'Done.'})
end


function pfits = fit_found_peak_core(data,opts)

    if isfield(opts,'fig_idx')
        fnum = 54714 + opts.fig_idx;
    else
        fnum = 54714;
    end
    f_cen = mean(data.peaks.freqs);
    sfigure(fnum);
    clf;
    subplot(2,1,1)
    plot(data.spec.freq-f_cen,data.spec.signal,'kx')
    hold on

    
    f_fit = data.spec.freq - f_cen;
    f_guesses = data.peaks.freqs-f_cen;
    h_guesses = data.peaks.vals;
    w_guesses = data.peaks.widths/2;
    gfun = @(p,x) 1-p(1)*exp(-((x-p(2))/p(3)).^2);
    lfun = @(p,x) 1-p(1)*p(3)./((x-p(2)).^2 + (p(3)).^2);
    g0 = [h_guesses,f_guesses,sqrt(w_guesses)];
    l0 = [0.5*w_guesses*(1-min(data.spec.signal)),f_guesses,0.5*w_guesses];
    
%     X_l  =data.spec.freq(2:end);
%     plot(data.spec.freq,gfun(g0,data.spec.freq),'rx')
    
    gfit = fitnlm(f_fit ,data.spec.signal,gfun,g0);
    lfit = fitnlm(f_fit ,data.spec.signal(1:end),lfun,l0);
    plot(f_fit ,gfun(gfit.Coefficients.Estimate,f_fit ),'bo')
    plot(f_fit ,lfun(lfit.Coefficients.Estimate,f_fit ),'ro')
    
%     plot(data.spec.freq,lfun(l0,data.spec.freq),'gx-')
    xlabel(strrep(data.data.class{1},'_',' '))
    pfits = [];
    
    
    legend('raw data','Gaussian fit','Lorentzian fit')
    
    subplot(2,1,2)
    plot(f_fit,gfun(gfit.Coefficients.Estimate,f_fit)-data.spec.signal,'bx')
    hold on
    plot(f_fit,lfun(lfit.Coefficients.Estimate,f_fit)-data.spec.signal,'rx')
    legend('Gaussian residuals','Lorentzian residuals')
    
    
end
