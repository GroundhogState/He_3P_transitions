% relative to 744396194.68
newpeaks{1}.data = [1.344,  2.673;
                   -12.816, -17.050;
                    2.127, 3.044];
newpeaks{2}.data = [3.621 0.270;
                    -1.336  -4.328;
                    3.857 0.914];
                
newpeaks{1}.ref = 744396208.36;
newpeaks{2}.ref = 744396208.36;
% Two-peak results STAGE 1:
%  - PEAK 1
% Strength 
% Centre   
% Width    
%  - PEAK 2
% Strength
% Centre   
% Width    
% Centres relative to 744396208.36
% DONE
% Two-peak results STAGE 2:
%  - PEAK 1
% Strength 
% Centre   
% Width    
%  - PEAK 2
% Strength 
% Centre  
% Width   
% Centres relative to 744396208.36
% DONE

% % FINAL PLOTS FOR D_2,3
lfun = @(p,x) p(1)*p(3)./((x-p(2)).^2 + (p(3)).^2);
ncat = numel(data.cat);
savedir = 'C:\Users\jaker\GoogleDrive\HEBEC\Thesis\fig\Spectroscopy';
cli_header(2,'Making plot...');
f1=stfig('Final plot mod');
clf;

data_colours = [226,75,42;
                45,225,245]/255;
theory_colours = 0.8*data_colours;
e_levels{1} = {'5^3D_3_2','5^3D_3_1','5^3D_2_2','5^3D_2_1'};
e_levels{2} = {'5^3D_3_2','5^3D_3_1','5^3D_2_2','5^3D_2_1'};
fmt_name = '5_3D_3';
f_pred = opts.const.f_table.g_2_3P_2.(sprintf('e_%s',fmt_name))/1e6;
% e_level

st_order = [1,3,2,4]; %(J mJ) = (32,31,22,21)
for cidx =1:ncat
    cdata = data.cat{cidx};
%     maxpt = max(cdata.spec.signal);
    normalizer = 1;
    B = opts.Bfield(cidx);
    subplot(2,1,cidx)
    hold on
    [~,p_cen] = max(cdata.pfits.lor_prms(:,3));
    f_offset = f_pred;
    
    all_X = linspace(min(X_data),max(X_data),500);
    X_data=cdata.spec.freq;
    signal = cdata.spec.signal;
    
    sat_mask = signal > 0.9;
    X_data = X_data(~sat_mask);
    signal = signal(~sat_mask);
    
    pastfreq = cdata.pfits.lor_prms(:,1);
    pastwidth = cdata.pfits.lor_prms(:,2);
%     pastheight = cdata.pfits.lor_prms(:,2);
    % s1, f1, w1
    guess1 = [1,pastfreq(1)-f_offset,pastwidth(1)];
    guess2 = [1,pastfreq(3)-f_offset,pastwidth(3)];
    newbits = newpeaks{cidx}.data';
    Delta = newpeaks{1}.ref - f_offset;
    
    sumspec = twopeak([newbits(1,:),newbits(2,:)],all_X-f_offset-Delta) + ...
                    twopeak([guess1,guess2],all_X-f_offset);
    prms = [newbits(1,:),newbits(2,:),guess1,guess2];
    fitprms = fitnlm(X_data-f_offset,signal/normalizer,@fourpeak,prms);
    prm_out = fitprms.Coefficients.Estimate;
    
%     plot(X_data-f_offset,fourpeak(prms,X_data-f_offset),'k')
    plot(all_X-f_offset,fourpeak(prm_out,all_X-f_offset),'k')
    plot(X_data-f_offset,signal/normalizer,'.','MarkerSize',15,'Color',data_colours(cidx,:))

    fprintf('Two-peak results STAGE %u:\n',cidx);
    fprintf('Centres: [%.2f,%.2f,%.2f,%.2f]\n',[prm_out(2),prm_out(5),prm_out(8),prm_out(11)]+f_offset);
    fprintf('relative to %.2f\n',f_offset);
    fprintf('Widths: [%.2f,%.2f,%.2f,%.2f]\n',2*[prm_out(3),prm_out(6),prm_out(9),prm_out(12)]);
    for pk = 1:4
        plot([prm_out(2+3*(pk-1))-1*prm_out(3+3*(pk-1)),prm_out(2+3*(pk-1))+1*prm_out(3+3*(pk-1))],[0.25,0.25],'r')
    end
    

        d_stage = zeros(2,5);
        d_stage(1,:) =[-347.8536;
                     -316.5440;
                     -309.8524;
                     -284.6747;
                      -54.7724];
        d_stage(2,:)=[ -330.4973
                     -310.7371
                     -308.0653
                     -289.6403
                      -39.7930];
        f_wrt = 744396511.00;
        delta = f_wrt-744396208.36;
        for thidx = 1:4
            % Compute theory error
            e_level = e_levels{cidx}{thidx};
            fmt_name = strrep(e_level,'^','_');
%             f_pred = opts.const.f_table.g_2_3P_2.(sprintf('e_%s',fmt_name))/1e6; %MHz
            B = opts.Bfield(cidx);
            % For the specified transition, compute the expected Zeeman shift
            g_state = '2_3P_2_2';
            g_level = g_state(1:6);
            %     zlines = z_state2state(B,g_state,strrep(e_state,'^','_'),const);
            g_e = lande_sg(e_level,opts.const);
            g_g = lande_sg(g_level,opts.const);
            m_g = str2double(strrep(strrep(g_state(end-1:end),'_',''),'n','-'));
            m_e = str2double(strrep(strrep(e_state(end-1:end),'_',''),'n','-'));
            del_mg = (g_e*m_e-g_g*m_g);
            E_shift = opts.const.mu*B*del_mg;
            f_shift = 1e-6*E_shift/opts.const.h; %MHz
            theory_error = 1e-6*opts.const.mu*del_mg*opts.Bfield_unc(cidx)/opts.const.h;
            f_shift = d_stage(cidx,thidx);
%             theory_error = .7; % indicative...
            X_thry = [-theory_error,-theory_error,theory_error,theory_error]+f_shift+delta;
            Y_thry = [0,1.2,1.2,0];
            thrybox=fill(X_thry,Y_thry,theory_colours(cidx,:),'FaceAlpha',0.3,'EdgeAlpha',0);
            uistack(thrybox,'bottom')
        end

        if cidx ==1
            set(gca,'FontSize',20)
            ylabel('Normalized  loss','FontSize',10)
            xticks([]);
        else
%             yticks([]);
            xlabel(sprintf('$\\nu$-%.2f (MHz)',f_offset))
        end
    set(gca,'FontSize',20)
    ylim([-0,1.2])
    
%     plot(Xfit-f_offset,all_fit/normalizer,'k:','LineWidth',2)
end
%%
savedir = 'C:\Users\jaker\GoogleDrive\HEBEC\Thesis\fig\Spectroscopy';
imname = sprintf('out_5_3D_23_pretty');
filename1 = fullfile(savedir,imname);
            saveas(f1,[filename1,'.fig']);
            saveas(f1,[filename1,'.png']);
            saveas(f1,[filename1,'.svg']);

cli_header(2,'Done');
% %
f2=stfig('Constraining levels');
clf;
hold on
Y1 = expdata.peaks(1,:)-f_23P2_53D1;
Y2 = expdata.peaks(2,:)-f_23P2_53D1;
X1 = expdata.field(1)*ones(size(Y1));
X2 = expdata.field(2)*ones(size(Y1));
Yerr = 2.2*ones(size(Y1));
Xerr = 0.3*ones(size(Y1));
th2=plot(plot_fields,(evals_opt)/1e6-f_23P(plot_fields)','k','LineWidth',1.5);
errorbar(X1,Y1,Yerr,Yerr,Xerr,Xerr,'.','MarkerSize',30,'MarkerFaceColor',data_colours(1,:));
errorbar(X2,Y2,Yerr,Yerr,Xerr,Xerr,'.','MarkerSize',30,'MarkerFaceColor',data_colours(2,:));
% th1=plot(plot_fields,sort(evals,2)/1e6-f_23P(plot_fields)','k:');
ylim([-370,-270])
set(gca,'FontSize',25)
% title('Constraining ','FontSize',35)
xlabel('Magnetic field strength (G)','FontSize',30)
ylabel(sprintf('$\\nu$ - %u (MHZ)',f_23P2_53D1),'FontSize',30)
imname = sprintf('out_constrain_splittings');
filename1 = fullfile(savedir,imname);
            saveas(f2,[filename1,'.fig']);
            saveas(f2,[filename1,'.png']);
            saveas(f2,[filename1,'.svg']);
cli_header(0,'Done!');
%% From the first attempt at fits
chidx = 2;
datachoice = data.cat{chidx};
freqs = datachoice.spec.freq;
ampl = datachoice.spec.signal;
normalizer = max(ampl);
f_ref = f_offset;

lfun = @(p,x) p(1)*p(3)./((x-p(2)).^2 + (p(3)).^2);
x_fit = linspace(-30,10,100);

pks_out = zeros(10,10,2);
% for cen_os1 = 1:10
%     for cen_os2 = 1:10

        % strength, centre, width
p1 = [1,10,3];
p2 = [1,10,3];

fit_out = fitnlm(freqs-f_ref,ampl/normalizer,@twopeak,[p1,p2]);
fprm = fit_out.Coefficients.Estimate;
pks_out(cen_os1,cen_os2,:) = fprm([2,5]);

stfig('khee');
clf;
subplot(2,1,1)
hold on
plot(freqs-f_ref,ampl/normalizer,'k.')
plot(x_fit,twopeak([p1,p2],x_fit))
plot(x_fit,twopeak(fprm,x_fit))
fprintf('Two-peak results STAGE %u:\n',chidx);
for pk = 1:2
    fprintf(' - PEAK %u\n',pk);
    fprintf('Strength %.3f\n',fprm(1+3*(pk-1)));
    fprintf('Centre   %.3f\n',fprm(2+3*(pk-1)));
    fprintf('Width    %.3f\n',fprm(3+3*(pk-1)));
end
xlabel(sprintf('Centres relative to %.2f\n',f_ref));
fprintf('Centres relative to %.2f\n',f_ref);
fprintf('DONE\n')
subplot(2,1,2)
hold on
histogram(squeeze(pks_out(:,:,1)))
histogram(squeeze(pks_out(:,:,2)));

function Y = twopeak(p,x)
    % s1, f1, w1
    f_cen1 = p(2);
    width1 = p(3);
    str1 = p(1);
    f_cen2 = p(5);
    width2 = p(6);
    str2 = p(4);
    C1 =  str1*width1./((x-f_cen1).^2 + (width1).^2);
    C2 =  str2*width2./((x-f_cen2).^2 + (width2).^2);
    Y = C1+C2;
end

function YY = fourpeak(p,x)
    YY = twopeak(p(1:6),x) + ...
                    twopeak(p(7:12),x);
end



