    % LV, Power, PD 
set = [0.1, 3.31, 0.0475;
        0.15, 5.14, 0.0645;
        0.2, 6.5, 0.0817;
        0.25, 8.9, 0.099];
nullset = [0, 0, 0.011];
delta = set-nullset;
mean(delta(:,2)./delta(:,3)) % conversion factor photodiode -> power (mW) ~ 95
% 1.1667 conversion from beam power -> peak power

% Line = Beam waist(cm) , Beam power (mW) ,  Peak intensity (W/m$^2$) , Exposure time (ms)];
triplet_S_1  =[ 4.1 ,  2.73 , 3.1 , 100 ];
triplet_D_1  =[ 4.1 ,  4.5 , 5.2 , 150];
triplet_D_23 =[ 4.1 ,  8.6 , 9.9 , 250  ];
singlet_D_2  =[ 1e-2 ,  10 , 6.3*10^3 , 100 ];

% Stark shift:

% h = 6.63e-34;
% c = 299792458;

G = 2*pi*1.6e6;
w0 = 2*pi*c.c/1083.331e-9;
f0 = c.c/1083.331e-9;
fscan = linspace(-5e6,5e6,500);
intens = 2*1e-8/(pi*1e-6);
alpha = @(Gamma,w,w_0) 6*pi*c.eps*c.c^3 * (Gamma/w_0^2)./(w_0^2 - w.^2 - 1j*Gamma*(w./w_0).^3);
df = @(I,a) -I*real(a)/(2*c.eps*c.c*c.h);

stfig('probe stark');
clf;
plot(fscan,real(alpha(G,fscan,w0)))


% % Doppler broadening

dopplerlimit = @(Gamma) c.h * Gamma / (2*1.308e-23);
dopplerlimit(1.6e6)
FWHM_broad = @(m,T) sqrt(8*T*1.308e-23*log(2)/(m*c.c^2));
FWHM_broad(6.64e-27,dopplerlimit(1.6e6))*7.441.6e14