function data = present_plots(data,opts)
    lfun = @(p,x) p(1)*p(3)./((x-p(2)).^2 + (p(3)).^2);
    cli_header({0,'Forming final plots...'});
    ncat = numel(data.cat);
    if ~iscell(opts.e_state)       
        num_pks = 1;
    else
        num_pks = numel(opts.e_state);
    end
    if isfield(opts,'fidx')
        fnum = 8888+opts.fidx;
    else
        fnum = 8888;
    end
    
   f1=sfigure(fnum);
   clf;
   [~,p_cen] = max(data.cat{1}.pfits.lor_prms(:,3));
%    f_offset = data.cat{1}.pfits.lor_prms(p_cen,1);
    f_offset = data.cat{1}.zeeman.f_pred(p_cen);
%    cmap = 
    for cidx =1:ncat
       
       cdata = data.cat{cidx};
       B = opts.Bfield(cidx);
%        fnum = 7777+cidx;
%        [~,p_cen] = max(cdata.pfits.lor_prms(:,3));
%        f_offset = cdata.pfits.lor_prms(p_cen,1);
       X_data=cdata.spec.freq;
       all_fit=zeros(1,150);
       if cidx ==1
           colour = [226,65,42]/255;
       else
           colour = [45,225,245]/255;
       end
      
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

            peak_df = cdata.zeeman.shift(pidx)-cdata.zeeman.shift(p_cen);
            f_cen = cdata.zeeman.corrected(p_cen);

            l_mdl = cdata.pfits.lorz{pidx};
            df = cdata.pfits.lor_prms(pidx,1) - cdata.pfits.lor_prms(p_cen,1);
            X=cdata.spec.freq;
            Xfit = linspace(min(X),max(X),150);
            Yfit = lfun(l_mdl.Coefficients.Estimate,Xfit-cdata.pfits.lor_prms(pidx,1));
            all_fit = all_fit+Yfit;
            [ysamp,yci]=predict(l_mdl,(Xfit-cdata.pfits.lor_prms(pidx,1))','Prediction','observation'); 

            freq_stat_err = cdata.zeeman.shift_unc(pidx) + cdata.pfits.lor_err(pidx,1);
            % Z shift unc + fit err 
            theory_error = 0.7+cdata.zeeman.shift_unc(pidx);
                % Drake's quoted error + zeeman shift error for plot
            % Plot theory values

%             if num_pks < 2,area(Xfit-f_offset,yci,'k:','LineWidth',1.);end
            if num_pks < 2
                x_boundary = [Xfit;Xfit]-f_offset;
                ar=area(x_boundary',yci,'FaceAlpha',0.16,'LineStyle','none');
                ar(1).FaceColor = [1 1 1];
                ar(2).FaceColor = 0.65*[1 1 1];
                plot([f_pred+f_shift,f_pred+f_shift]-f_offset,[-1e4,0],'Color',[0.8500 0.3250 0.0980])
                plot([f_pred+f_shift,f_pred+f_shift]-f_offset+theory_error,-1e4*[.6,.4],'Color',[0.8500 0.3250 0.0980])
                plot([f_pred+f_shift,f_pred+f_shift]-f_offset-theory_error,-1e4*[.6,.4],'Color',[0.8500 0.3250 0.0980])
                plot([f_pred+f_shift-f_offset-theory_error,f_pred+f_shift-f_offset+theory_error],-1e4*[.5,.5],'Color',[0.8500 0.3250 0.0980])
                ms = 12;
                supt = sprintf('%s%s','2^3P_2\rightarrow',e_level);
            else
                ms = 15;
                supt = sprintf('%s5^3D_2,%s','2^3P_2\rightarrow',e_level);
            end
%             plot(Xfit-f_offset,Yfit,'kx')
        end
         plot(X_data-f_offset,cdata.spec.signal,'.','MarkerSize',ms,'Color',colour)
         plot(Xfit-f_offset,all_fit,'k')
    end
           
            ylabel('\Delta')
            ylim(1e4*[-1,7])
%             xlabel('\nu')
            xlabel(sprintf('\\nu-%.2f (MHz)',f_offset))
%             legend({'Data','Theoretical value'},'location','NorthWest')
            set(gca,'FontSize',15,'FontName',opts.font)
            set(gcf,'color','w');
            box off
%             XAxisLocation = 'bottom'
            st=suptitle(supt);
            st.FontName = 'Times';
            imname = sprintf('%s_plot_combo',e_level);
            filename1 = fullfile(opts.out_dir,imname);
            saveas(f1,strcat(filename1,".fig"));
            saveas(f1,strcat(filename1,".svg"));
            saveas(f1,strcat(filename1,".png"))
            
end