function fits = fit_found_peaks(data,opts)

    pks = data.tr.peaks;
    num_pks = length(pks.vals);
    header({0,'Fitting %u detected peaks',num_pks})

    % Probably need to mask up the data to help fit the peaks... OR fit the combination thereof?!
    Xdata = data.spec.freq;
    Ydata = 1-data.spec.signal;    
    Xmdl = linspace(min(Xdata),max(Xdata),250);
    
    gaus_mdlfun = @(p,x) p(1)*exp(-0.5*((x-p(2))/p(3)).^2);
    lorz_mdlfun = @(b,x) b(1).*b(3)./((x-b(2)).^2+b(3).^2);
    
    sfigure(4752);
    clf;
    plot(Xdata,Ydata,'ko')
    hold on
    
    for pidx = 1:num_pks
        pk = struct_mask(pks,pidx);
        g0 = [pk.vals,pk.freqs,pk.widths];
        l0 = [pi*pk.vals,pk.freqs,0.5*pk.widths];
        
        %% Gaussian weights that fall off like the peak
        w = exp( -((Xdata - pk.freqs)./(0.5*pk.widths)).^2)+eps;
        
        plot(Xmdl,gaus_mdlfun(g0,Xmdl),'r.')
        plot(Xmdl,lorz_mdlfun(l0,Xmdl),'k.')
        plot(Xdata,w,'g-')
        
        gfit = fitnlm(Xdata,Ydata,gaus_mdlfun,g0,'Weight',w);
        lfit = fitnlm(Xdata,Ydata,lorz_mdlfun,l0,'Weight',w); 
        plot(Xmdl,gaus_mdlfun(gfit.Coefficients.Estimate,Xmdl),'r-')
        plot(Xmdl,lorz_mdlfun(lfit.Coefficients.Estimate,Xmdl),'k-')
        
    end
    
    legend('Data','Gaussian guess,','Lorentzian guess','Weighting fn','Gaussian fit','Lorentzian fit')

    fits = [];

    header({1,'Done.'})
    
end