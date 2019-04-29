function z_correct = zeeman_correction(data,opts)
    header({0,'Correcting for Zeeman shift'})
    numcat = numel(data.cat);
    f_shifted = zeros(numcat,1);
    for cidx=1:numcat
        B = opts.Bfield(cidx);
        % For the specified transition, compute the expected Zeeman shift
        g_state = '2_3P_2_2';
        g_level = g_state(1:6);
        e_state_format=strrep(e_state,'^','_');
        e_level = e_state_format(1:6);
    %     zlines = z_state2state(B,g_state,strrep(e_state,'^','_'),const);
        g_e = lande_sg(e_level,opts.const);
        g_g = lande_sg(g_level,opts.const);
        m_g = str2num(strrep(strrep(g_state(end-1:end),'_',''),'n','-'));
        m_e = str2num(strrep(strrep(e_state(end-1:end),'_',''),'n','-'));
        del_mg = (g_e*m_e-g_g*m_g);
        E_shift = opts.const.mu*B*del_mg;
        f_shift = 1e-6*E_shift/opts.const.h; %MHz
        f_shifted(cidx) = data.cat{cidx}.pfits.lor_prms(:,1) - f_shift;
        fprintf('Z-corrected centre freq %u:       %.2f(%.2f) MHz\n',cidx,f_shifted(cidx),data.cat{cidx}.pfits.lor_err(1))
    end 
    z_correct.f_shifted = f_shifted;
end
