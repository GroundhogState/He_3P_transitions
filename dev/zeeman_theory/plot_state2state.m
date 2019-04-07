function plot_state2state(lines,g_state,e_state,const,opt)
            g_level = g_state(1:6);
            e_level = e_state(1:6);
            e_term = e_level(1:4);
%             g_man = g_term(1:3);
            e_man = e_term(1:3);
            m_g = strrep(g_state(end-1:end),'_','');
            m_e = strrep(e_state(end-1:end),'_','');
            mg_num = str2num(strrep(m_g,'n','-'));
            me_num = str2num(strrep(m_e,'n','-'));
            J_e = e_level(end);
            level_pair = lines.(e_term(4)).(['g_',g_level]).(['e_',e_term,'_',num2str(J_e)]);
            state_pair = level_pair.(['mg_',m_g]).(['me_',m_e]);
            exclude_transition = false;
            switch me_num-mg_num
                case -1
                    lstyle = 'v';
                case 0
                    lstyle = 'o';
                case 1
                    lstyle = '^';
                otherwise
                    exclude_transition = true; %if not allowed by dipole rules
            end
            if ~exclude_transition
                p=plot(state_pair.B,state_pair.f,lstyle);
                if isfield(opt,'plot_colour')
                    p.Color=opt.plot_colour;
                end
                hold on
            end
end