function manifold_lines = z_level2manifold(B,g_level,e_manifold,const)

    L_set = 0:2; %hard coding for now
    for lidx = 1:length(L_set)
%         L_e = L_set(lidx);
        e_term = [e_manifold,const.terms{lidx}]
        manifold_lines.(const.terms{lidx}) = z_level2term(B,g_level,e_term,const)
    end
    
end