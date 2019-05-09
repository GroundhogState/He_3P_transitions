function data = zeeman_correction(data,opts)
    if isfield(opts,'e_state')
        cli_header({0,'Correcting for Zeeman shift'})
        
%         e_state_format=strrep(opts.e_state,'^','_');
        
        if ~iscell(opts.e_state)       
            num_pks = 1;
        else
            num_pks = numel(opts.e_state);
        end

        
        
        
        numcat = numel(data.cat);
        
%         f_shifted = zeros(numcat,1);
        % This is kinda hacky just because we're lucky enough that the
        % first peak is the one we wanna fit
        for cidx=1:numcat
            for npk=1:num_pks
                if num_pks == 1
                    e_state = opts.e_state;
                    e_level = e_state(1:6);
                    
                else
                    e_state = opts.e_state{npk};
                    e_level = e_state(1:6);
                end
                fmt_name = strrep(e_level,'^','_');
                f_pred = opts.const.f_table.g_2_3P_2.(sprintf('e_%s',fmt_name))/1e6; %MHz
                B = opts.Bfield(cidx);
                % For the specified transition, compute the expected Zeeman shift
                g_state = '2_3P_2_2';
                g_level = g_state(1:6);

            %     zlines = z_state2state(B,g_state,strrep(e_state,'^','_'),const);
                g_e = lande_sg(e_level,opts.const);
                g_g = lande_sg(g_level,opts.const);
                m_g = str2num(strrep(strrep(g_state(end-1:end),'_',''),'n','-'));
                m_e = str2num(strrep(strrep(e_state(end-1:end),'_',''),'n','-'));
                del_mg = (g_e*m_e-g_g*m_g);
                E_shift = opts.const.mu*B*del_mg;
                f_shift = 1e-6*E_shift/opts.const.h; %MHz
                f_shift_unc = 1e-6*opts.const.mu*del_mg*opts.Bfield_unc(cidx)/opts.const.h;
                f_shifted = data.cat{cidx}.pfits.lor_prms(npk,1) - f_shift;
                fprintf('Peak %u shifted %.1f MHz to:       %.2f(%.2f) MHz\n',...
                    npk,-f_shift,f_shifted,data.cat{cidx}.pfits.lor_err(npk)+f_shift_unc)
                fprintf('       theory val:                 %.2f(diff %.2f)MHz,  \n',...
                    f_pred,f_shifted-f_pred)
                data.cat{cidx}.zeeman.f_pred(npk) = f_pred;
                data.cat{cidx}.zeeman.shift(npk) = f_shift;
                data.cat{cidx}.zeeman.shift_unc(npk) = f_shift_unc;
                data.cat{cidx}.zeeman.corrected(npk)= f_shifted;
                data.cat{cidx}.zeeman.g_factors(npk,:) = [g_g,g_e];
                data.cat{cidx}.zeeman.stat_unc(npk) = data.cat{cidx}.pfits.lor_err(npk);
            end
        end
        
    else
        cli_header({1,'Zeeman correction bypassed'})
        z_correct = [];
    end
end
