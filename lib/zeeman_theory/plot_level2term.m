

function plot_level2term(lines,g_level,e_term,const,opt)
    J_g = str2num(g_level(end));
    if ~isfield(opt,'mg_range')
        mg_range = -J_g:J_g; % options to restric
        mg_range = 2;
    end
    for i=1:length(mg_range)
%         cmap = colormap(viridis(length(mg_range)+1));
%         opt.plot_colour = cmap(i,:);
        mg = mg_range(i);
        g_state = [g_level,'_',strrep(num2str(mg),'-','n')];
        plot_state2term(lines,g_state,e_term,const,opt)
    end
end
