function g = lande_sg(level,const)
    % Computes the g-factor from spectroscopic notation
    N = str2num(level(1));
    S = 0.5*(str2num(level(3))-1);
    [~ ,i] = find(strcmp(level(4),const.terms));
    L = i-1;
    J=str2num(level(6));
    sprintf('S=%u, L=%u, J=%u',S,L,J); % debug
    g = lande_g(S,L,J);
end