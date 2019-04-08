function plot_state2term(lines,g_state,e_term,const,opt)
    S_e = 0.5*(str2num(e_term(3))-1);
    [~ ,i] = find(strcmp(e_term(4),const.terms));
    L_e = i-1;
    J_e_range = abs(L_e - S_e):L_e + S_e;%for example, selection rules come later
    for J_idx = 1:length(J_e_range)
        cmap = colormap(viridis(length(J_e_range)+1));
        opt.plot_colour = cmap(J_idx,:);
        J_e = J_e_range(J_idx);
        e_level = [e_term,'_',num2str(J_e)];
        plot_state2level(lines,g_state,e_level,const,opt)
    end
end