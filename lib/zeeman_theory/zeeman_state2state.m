function state_lines = zeeman_state2state(B,g_state,e_state,const)
    % Calculate g-factors
    g_level = g_state(1:6);
    e_level = e_state(1:6);
%     g_e = lande_sg(g_level,const);
%     g_g = lande_sg(e_level,const);
    g_e = lande_sg(e_level,const);
    g_g = lande_sg(g_level,const);
    m_g = str2num(strrep(strrep(g_state(end-1:end),'_',''),'n','-'));
    m_e = str2num(strrep(strrep(e_state(end-1:end),'_',''),'n','-'));
    % Compute the absolute & differential energies & frequencies
    del_mg = g_e*m_e-g_g*m_g;
    f0 = const.f_table.(['g_',g_level]).(['e_',e_level]);
    E0 = const.h*f0;
    E_shift = const.mu*B*del_mg;
    E_z = E0+E_shift;
    state_lines.B = B;
    state_lines.E = E_z;
    state_lines.f = E_z/const.h;
    state_lines.dE = E_shift;
    state_lines.df = E_shift./const.h;
    state_lines.wl = (const.h*const.c)./E_z;
    state_lines.d_wl = (const.h*const.c)./E_shift;
end