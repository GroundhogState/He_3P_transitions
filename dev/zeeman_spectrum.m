% Let's start with the big boy; 2^3P_2 to 5^3D_3
disp('=== Zeeman splittings ===')
clear all;
% Physical constants
const.mu = 9.27e-28; %J/G
const.h = 6.63e-34;
const.hbar = const.h/(2*pi);
const.f_mu = const.mu/const.h;
const.w_mu = const.mu/const.hbar;
const.c = 299792458;
% Notation & lookup
const.terms = {'S','P','D','F','G'};
% Use underscores for valid field names
% Read N_(2S+1)L_J 
const.f_table.g_2_3P_2.e_5_3D_3 = 744.39620968e12; %g, e, f0 (Hz)
const.f_table.g_2_3P_2.e_5_3D_2 = 744.39622889e12;
const.f_table.g_2_3P_2.e_5_3D_1 = 744.39651246e12; 
const.f_table.g_2_3P_2.e_5_1D_2 = 744.43034335e12;
% Note that J ranges from |L-S| to L+S
% A combination of S and L is called a term
% S, L, J specify a level
% S, L, J, and m_J specify a state
% Parameters
B = 1;
g = '2_3P_2';
e = '5_3D_3';

S = zeeman_splitting(B,g,e,const);
freqs = S(:,3);

lower_level = '2_3P_2';
upper_term = '5_3D';
prm.Brange = 0:25;
prm.Bshow = 10;
prm.f_plot_offset = 7.44396e8;

plot_zeeman_fan(lower_level,upper_term,prm,const)

prm.f_plot_offset = 0;
plot_zeeman_spectrum(lower_level,upper_term,prm,const)

prm.


disp('  = Done =  ')


%% functions
function plot_zeeman_spectrum(lower_level,upper_term,prm,const)
    figure(2)
    clf;
    S_e = 0.5*(str2num(upper_term(3))-1);
    [~ ,i] = find(strcmp(upper_term(4),const.terms));
    L_e = i-1;
    J_vals = abs(L_e-S_e):L_e+S_e;
    for m=J_vals
        upper_level = [upper_term,'_',num2str(m)];
        spectrum = zeeman_splitting(prm.Bshow ,lower_level,upper_level,const);
        mjs = spectrum(:,1:2);
        d_mj = mjs(:,2)-mjs(:,1);
        sigma_plus = d_mj==1;
        pi_ = d_mj ==0;
        sigma_minus = d_mj == -1;
        freqs = spectrum(:,3)/1e6;

        [X_plus,Y_plus]=vlines(freqs(sigma_plus)-prm.f_plot_offset,0,-0.2);
        plot(X_plus,Y_plus,'r')
        hold on
        [X_minus,Y_minus]=vlines(freqs(sigma_minus)-prm.f_plot_offset,0,-0.2);
        plot(X_minus,Y_minus,'b')
        [X_pi,Y_pi]=vlines(freqs(pi_)-prm.f_plot_offset,0,-0.2);
        plot(X_pi,Y_pi,'k')
    end
    lower_label = lower_level;
    lower_label(2) = '^';
    ylim([-0.2,1.1])
    title(sprintf('Peak positions of all %s-%s transitions at %.1f Gauss',lower_label,strrep(upper_term,'_','^'),prm.Bshow))
end

function [X,Y]= vlines(locs,ymin,ymax)
% accepts COLUMN vector of X 
    X=(locs.*[1,1])';
    Y=[ymin*ones(size(locs)),ymax*ones(size(locs))]';
end

function plot_zeeman_fan(lower_level,upper_term,prm,const)
    figure(1)
    clf
    num_B = length(prm.Brange);
    S_e = 0.5*(str2num(upper_term(3))-1);
    [~ ,i] = find(strcmp(upper_term(4),const.terms));
    L_e = i-1;
    J_vals = abs(L_e-S_e):L_e+S_e;
    for m=J_vals
        upper_level = [upper_term,'_',num2str(m)];
        zlines = [];
        for n=1:num_B
           zlines(n,:,:) = zeeman_splitting(prm.Brange(n),lower_level,upper_level,const); 
        end
        all_lines = zlines(:,:,3);
        % Show the spectrum
        plot(prm.Brange,all_lines/1e6-prm.f_plot_offset,'k:')
        hold on
    end
    lower_label = lower_level;
    lower_label(2) = '^';
    title(sprintf('Zeeman splitting of %s-%s transitions',lower_label,strrep(upper_term,'_','^')))
    xlabel('Magnetic field (Gauss)')
    ylabel(sprintf('Frequency -%.2f (MHz)',prm.f_plot_offset))
end





function spectrum = zeeman_splitting(B,g,e,const)
    % Returns a table of format [m_g,m_e,f]

    % Math
    % Calculate g-factors
    g_e = lande_sg(g,const);
    g_g = lande_sg(e,const);
    % Populate the lists of projection quantum numbers
    m_g = -str2num(g(end)):str2num(g(end));
    m_e = -str2num(e(end)):str2num(e(end));
    % Table of m_j changes
    del_mj = m_e'-m_g;
    % Changes in J_z
    del_mg = g_e*m_e'-g_g*m_g; %column = init m
    dipole_mask = abs(del_mj) <=1;
    [i,j,~] = find(dipole_mask);

    % Physics
    % Calculations
    f0 = const.f_table.(['g_',g]).(['e_',e]);% Theory value; will eventually look it up
    E0 = const.h*f0;
    % Energy shfits
    delta_E = -const.mu*B*del_mg;
    dipole_shift = delta_E.*dipole_mask;
    dipole_transitions = E0+dipole_shift;
    % Frequency shifts
    dipole_f_shift = dipole_shift/const.h;
    dipole_f_transitions = dipole_transitions/const.h;
    [~,~,dipole_shift_list] = find(dipole_f_transitions.*dipole_mask);
    spectrum = [m_g(j)',m_e(i)',dipole_shift_list];
%     spectrum = dipole_f_shift(dipole_mask);
end


function g = lande_g(S,L,J)
    g = 3/2 + (S*(S+1)-L*L(+1))/(2*J*(J+1));
end

function g = lande_sg(state,const)
    % Computes the g-factor from spectroscopic notation
    N = str2num(state(1));
    S = 0.5*(str2num(state(3))-1);
    [~ ,i] = find(strcmp(state(4),const.terms));
    L = i-1;
    J=str2num(state(end));
    sprintf('S=%u, L=%u, J=%u',S,L,J); % debug
    g = lande_g(S,L,J);
end


