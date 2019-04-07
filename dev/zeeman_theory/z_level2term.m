function term_lines = z_level2term(B,g_level,e_term,const)
% Returns all transitions between a specified initial state and a target
% term
    S_e = 0.5*(str2num(e_term(3))-1);
    [~ ,i] = find(strcmp(e_term(4),const.terms));
    L_e = i-1;
    J_vals = abs(L_e-S_e):L_e+S_e;
    % Loop over the levels
    for J = J_vals
        e_level = [e_term,'_',num2str(J)];
        term_lines.(['g_',g_level]).(['e_',e_level]) = zeeman_level2level(B,g_level,e_level,const);
    end
end