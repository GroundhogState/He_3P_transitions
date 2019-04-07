function level_lines = zeeman_level2level(B,g_level,e_level,const)
    % Returns all transitions between two levels
    % accepts a lower and upper term symbol; produces lines labeled by states 
    % Returns a table of format [m_g,m_e,f] of transitions between states
    
    % Calculate g-factors
    g_e = lande_sg(g_level,const);
    g_g = lande_sg(e_level,const);
    % Populate the lists of projection quantum numbers
    m_g_set = -str2num(g_level(end)):str2num(g_level(end));
    m_e_set = -str2num(e_level(end)):str2num(e_level(end));
    % Loop over the states
    for m_g=m_g_set
        for m_e = m_e_set
            g_state = [g_level,'_',strrep(num2str(m_g),'-','n')];
            e_state = [e_level,'_',strrep(num2str(m_e),'-','n')];
            
            g_level = g_state(1:6);
            e_level = e_state(1:6);
            if isfield(const.f_table,['g_',g_level])
                if isfield(const.f_table.(['g_',g_level]),['e_',e_level])
                    state_lines = zeeman_state2state(B,g_state,e_state,const);
                    level_lines.(['mg_',strrep(num2str(m_g),'-','n')]).(['me_',strrep(num2str(m_e),'-','n')]) = state_lines;
                end
            end
        end
    end

end