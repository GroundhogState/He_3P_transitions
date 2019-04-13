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
    sfigure(fnum);
    clf;
    plot(data.spec.freq,data.spec.signal,'k.')
    hold on

    
    f_guesses = data.peaks.freqs;
    h_guesses = data.peaks.vals;
    w_guesses = data.peaks.widths;
    npeaks = length(h_guesses);
    gfun = @(p,x) p(1)*exp(-((x-p(2))/p(3)).^2);
    b0 = zeros(npeaks,3);
    for nn=1:npeaks
        b0(nn,:) = [h_guesses(nn),f_guesses(nn),w_guesses(nn)];
    end
    S = spectral_func(b0,data.spec.freq,gfun);   
    plot(data.spec.freq,1-S,'o')
%     s_fun = @(P,X) spectral_func(P,X,gfun);
%     S_fit = fitnlm(data.spec.freq,data.spec.signal,@(P,X) spectral_func(P,X,gfun),b0);
%     plot(data.spec.freq,1-S_fit,'x')
    xlabel(data.data.class{1})
    pfits = [];
        
    
end

function S = spectral_func(P,X,fun)
    % Returns a sum of peaks with parameters specified by the ROWS P over domain X
    % Need to return a COLUMN vector to match Y
    if size(X,2) > size(X,1)
        X = X';
    end
    S = zeros(size(X));
    npeaks = size(P,1);
    for n=1:npeaks
        S = S + fun(P(n,:),X);
    end
end



