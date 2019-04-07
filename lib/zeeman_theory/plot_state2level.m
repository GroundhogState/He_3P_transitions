function plot_state2level(lines,g_state,e_level,const,opt)
    J_g = str2num(g_state(6));
    J_e = str2num(e_level(6));
    % declare which d_mj are allowed
        for me_num = max(-J_e,J_g-1):min(J_e,J_g+1) % Dipole selection rules for now
            m_e = strrep(num2str(me_num),'-','n');
            e_state = [e_level,'_',m_e];
            plot_state2state(lines,g_state,e_state,const,opt)
        end
end