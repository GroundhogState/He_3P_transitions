function data = present_plots(data,opts)
    lfun = @(p,x) p(1)*p(3)./((x-p(2)).^2 + (p(3)).^2);
    cli_header({0,'Forming final plots...'});
    ncat = numel(data.cat);
    if ~iscell(opts.e_state)       
        num_pks = 1;
    else
        num_pks = numel(opts.e_state);
    end
   f1=sfigure(7777);
   clf;
    for cidx =1:ncat
       
       cdata = data.cat{cidx};
       B = opts.Bfield(cidx);
%        fnum = 7777+cidx;
       subplot(2,1,cidx)
       [~,p_cen] = max(cdata.pfits.lor_prms(:,3));
       f_offset = cdata.pfits.lor_prms(p_cen,1);
       X_data=cdata.spec.freq;
       all_fit=zeros(1,100);
       plot(X_data-f_offset,cdata.spec.signal,'.','MarkerSize',15,'Color',[0.8350 0.1780 0.1840])
       hold on
        for pidx=1:num_pks
            if num_pks ==1
                e_state = opts.e_state;
            else
                e_state = opts.e_state{pidx};
            end
            e_level = e_state(1:6);
            fmt_name = strrep(e_level,'^','_');
            f_pred = opts.const.f_table.g_2_3P_2.(sprintf('e_%s',fmt_name))/1e6; %MHz
            f_shift = cdata.zeeman.shift(pidx);
            cli_header({1,'\nRESULTS for %s transition:\n',e_state})
            fprintf('Predicted vacuum freq %.3f MHz\n',f_pred)

            peak_df = cdata.zeeman.shift(pidx)-cdata.zeeman.shift(p_cen);
            f_cen = cdata.zeeman.corrected(p_cen);

            l_mdl = cdata.pfits.lorz{pidx};
            df = cdata.pfits.lor_prms(pidx,1) - cdata.pfits.lor_prms(p_cen,1);
            X=cdata.spec.freq;
            Xfit = linspace(min(X),max(X),100);
            Yfit = lfun(l_mdl.Coefficients.Estimate,Xfit-cdata.pfits.lor_prms(pidx,1));
            all_fit = all_fit+Yfit;
            [ysamp,yci]=predict(l_mdl,(Xfit-cdata.pfits.lor_prms(pidx,1))','Prediction','observation'); 

            freq_stat_err = cdata.zeeman.shift_unc(pidx) + cdata.pfits.lor_err(pidx,1);
            % Z shift unc + fit err 
            theory_error = 0.7+cdata.zeeman.shift_unc(pidx);
                % Drake's quoted error + zeeman shift error for plot

            fprintf('Peak %u fitted value:   %.2f(%.2f) MHz\n',pidx,cdata.zeeman.corrected(pidx),freq_stat_err)
            fprintf('         theory diff:   %.2f MHz\n',cdata.zeeman.corrected(pidx)-f_pred)

            plot([f_pred+f_shift,f_pred+f_shift]-f_offset,[-1e4,0],'Color',[0.8500 0.3250 0.0980])
            plot([f_pred+f_shift,f_pred+f_shift]-f_offset+theory_error,-1e4*[.6,.4],'Color',[0.8500 0.3250 0.0980])
            plot([f_pred+f_shift,f_pred+f_shift]-f_offset-theory_error,-1e4*[.6,.4],'Color',[0.8500 0.3250 0.0980])
            plot([f_pred+f_shift-f_offset-theory_error,f_pred+f_shift-f_offset+theory_error],-1e4*[.5,.5],'Color',[0.8500 0.3250 0.0980])
            if num_pks < 2,plot(Xfit-f_offset,yci,'k:','LineWidth',1.);end
%             plot(Xfit-f_offset,Yfit,'b')



        end
            plot(Xfit-f_offset,all_fit,'b')
            ylabel('Atoms lost')
            xlabel(sprintf('Frequency-%.2f (MHz)',f_offset))
            title(sprintf('|B|=%.2f',B))
            legend({'Data','Theoretical value'},'location','NorthWest')
            set(gca,'FontSize',10,'FontName','cmr10')
            set(gcf,'color','w');



    end
            suptitle(sprintf('%s transition',e_level))
            imname = sprintf('%s_plot_pretty_%u',e_level,cidx);
            filename1 = fullfile(opts.out_dir,imname);
            saveas(f1,[filename1,'.fig']);
            saveas(f1,[filename1,'.png'])
end