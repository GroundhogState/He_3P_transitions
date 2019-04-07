function plot_state2manifold(lines,g_state,e_manifold,const,opt)
    L_set = 1:3; %hard coding for now
    for lidx = 1:length(L_set)
        L_e = L_set(lidx);
        e_term = [e_manifold,const.terms{L_e}];
        plot_state2term(lines,g_state,e_term,const,opt)
    end
end