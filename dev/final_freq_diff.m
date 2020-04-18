
 %%
w_5_3S_2__J1 = 4121.973334e-10;
w_5_3D_2__J3 = 4027.323811e-10;
w_5_3D_2__J2 = 4027.323707e-10;
w_5_3D_2__J1 = 4027.322173e-10;
w_5_1D_1__J2 = 4027.151538e-10;% + 2291.177 MHz

f_5_3S_2__J1 = opts.const.c/w_5_3S_2__J1;
f_5_3D_2__J3 = opts.const.c/w_5_3D_2__J3;
f_5_3D_2__J2 = opts.const.c/w_5_3D_2__J2;
f_5_3D_2__J1 = opts.const.c/w_5_3D_2__J1;
f_5_1D_1__J2 = opts.const.c/w_5_1D_1__J2+2291.177e6; %from Zheng, as Drake only lists lines from the P J=1
cli_header('Final freqs (drake)');
fprintf('3S_1 %.2f\n',f_5_3S_2__J1/1e6)
fprintf('3D_3 %.2f\n',f_5_3D_2__J3/1e6)
fprintf('3D_2 %.2f\n',f_5_3D_2__J2/1e6)
fprintf('3D_1 %.2f\n',f_5_3D_2__J1/1e6)
fprintf('1D_2 %.2f\n',f_5_1D_1__J2/1e6)

cli_header(1,'Widths from fits');
wvals = [3.21,0.4,9.21;
        5.79,0.62,16.4;
        4.18,0.5,16.4;
        4.04,0.12,16.4;
        3.12,0.13,13.9];
    
 mean(wvals(:,1)./wvals(:,3));
 wvals(:,2)./wvals(:,1);
 
 % 
 f_meas = [727303249,744396496,744396217,744396201,744430347];
 f_pred = [f_5_3S_2__J1,f_5_3D_2__J1,f_5_3D_2__J2,f_5_3D_2__J3,f_5_1D_1__J2]/1e6;
 df = f_meas-f_pred