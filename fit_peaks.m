function fits = fit_peaks(data,opts)

    [peak_guess,peak_idx]= min(data.tr.freq_stats.sig_cal);
    peak_freq = data.tr.freq_stats.freq(peak_idx);
    
    
    
    freqs_centred_fit = data.tr.sync.msr.probe_set-peak_freq;
    var_guess = 15;
    % Fit a gaussian
    gaus_mdlfun = @(p,x) 1 - p(1)*exp(-0.5*((x-p(2))/p(3)).^2);
    beta0 = [1-peak_guess,0,var_guess];
    fit_gaus = fitnlm(freqs_centred_fit,data.tr.sync.msr.N_atoms./data.tr.sync.msr.calib',gaus_mdlfun,beta0);

    % Fit Lorentzian
    var_guess = 15;
    lor_mdlfun = @(b,x) b(1)./((x-b(2)).^2+b(3));
    lor_mdlfun = @(b,x) b(1)./((x-b(2)).^2+b(3))+1;
    beta0 = [-peak_guess,0,var_guess];
    fit_lor = fitnlm(freqs_centred_fit,data.tr.sync.msr.N_atoms./data.tr.sync.msr.calib',lor_mdlfun,beta0);
    lor_prms = fit_lor.Coefficients.Estimate;
    % Voigt profile?
    
    
    
    %extract useful parameters
    lorz_cen = fit_lor.Coefficients.Estimate(2)+peak_freq;
    lorz_unc = fit_lor.Coefficients.SE(2);
    lorz_width = 2*sqrt(abs(fit_lor.Coefficients.Estimate(3)));
    lorz_width_unc = 2*fit_lor.Coefficients.Estimate(3);
    lorz_wav = (299792458/lorz_cen)*1e3;
    lorz_cen_Hz = lorz_cen*1e6;
    lorz_unc_Hz = lorz_unc*1e6;
    lorz_wav_unc = lorz_unc_Hz*299792458/((lorz_cen_Hz)^2) * 1e9; % convert back to nm
    

    gaus_cen = fit_gaus.Coefficients.Estimate(2)+peak_freq;%+opts_tr.aom_freq;
    gaus_unc = fit_gaus.Coefficients.SE(2);
    gaus_width = fit_gaus.Coefficients.Estimate(3);
    gaus_width_unc = fit_gaus.Coefficients.SE(3);
    gaus_cen_Hz = lorz_cen*1e6;
    gaus_unc_Hz = lorz_unc*1e6;
    gaus_wav = (299792458/gaus_cen)*1e3;
    gaus_wav_unc = gaus_unc_Hz*299792458/((gaus_cen_Hz)^2) * 1e9; % convert back to nm
    
    % 
    
    
    gaus_prms = fit_gaus.Coefficients.Estimate;
    gaus_cen_prms = gaus_prms;
    gaus_cen_prms(2) = gaus_cen_prms(2) - fit_lor.Coefficients.Estimate(2);
    lor_cen_prms = lor_prms;
    lor_cen_prms(2) = 0;
    
    
    fits.freqs_centred = data.tr.sync.msr.probe_set-lorz_cen;
    fits.lor_mdlfun = lor_mdlfun;
    fits.lor_cen_prms = lor_cen_prms;
    fits.gaus_mdlfun = gaus_mdlfun;
    fits.gaus_cen_prms = gaus_cen_prms;
    fits.gaus_fit = [gaus_cen,gaus_unc,gaus_width,gaus_width_unc,gaus_wav,gaus_wav_unc];
    fits.lorz_fit = [lorz_cen,lorz_unc,lorz_width,lorz_width_unc,lorz_wav,lorz_wav_unc];
    
end

