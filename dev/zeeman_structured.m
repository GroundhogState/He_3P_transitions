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
B = 1:10;
g_level = '2_3P_2';
e_term = '5_3D';


lines = z_level2term(B,g_level,e_term,const)

disp('  = Done =  ')




%% functions
function term_lines = z_level2term(B,g_level,e_term,const)
% Returns all transitions between a specified initial state and a target
% term
    S_e = 0.5*(str2num(e_term(3))-1);
    [~ ,i] = find(strcmp(e_term(4),const.terms));
    L_e = i-1;
    J_vals = abs(L_e-S_e):L_e+S_e;
    for J = J_vals
        e_level = [e_term,'_',num2str(J)];
        term_lines.(['g_',g_level]).(['e_',e_level]) = zeeman_level2level(B,g_level,e_level,const);
    end
end

function transition_lines = zeeman_level2level(B,g_level,e_level,const)
% Returns all transitions between two levels
    % accepts a lower and upper term symbol; produces lines labeled by states 
    % Returns a table of format [m_g,m_e,f] of transitions between states
    
    % Math
    % Calculate g-factors
    g_e = lande_sg(g_level,const);
    g_g = lande_sg(e_level,const);
    % Populate the lists of projection quantum numbers
    m_g_set = -str2num(g_level(end)):str2num(g_level(end));
    m_e_set = -str2num(e_level(end)):str2num(e_level(end));
    % Generate the lookup table
    for m_g=m_g_set
        for m_e = m_e_set
            del_mg = g_e*m_e-g_g-m_g;
            f0 = const.f_table.(['g_',g_level]).(['e_',e_level]);
            E0 = const.h*f0;
            E_shift = -const.mu*B*del_mg;
            E_z = E0+E_shift;
            transition_lines.(['mg_',strrep(num2str(m_g),'-','m')]).(['me_',strrep(num2str(m_e),'-','n')]).E = E_z;
            transition_lines.(['mg_',strrep(num2str(m_g),'-','m')]).(['me_',strrep(num2str(m_e),'-','n')]).F = E_z/const.h;
            transition_lines.(['mg_',strrep(num2str(m_g),'-','m')]).(['me_',strrep(num2str(m_e),'-','n')]).dE = E_shift;
            transition_lines.(['mg_',strrep(num2str(m_g),'-','m')]).(['me_',strrep(num2str(m_e),'-','n')]).dF = E_shift/const.h;
        end
    end

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


