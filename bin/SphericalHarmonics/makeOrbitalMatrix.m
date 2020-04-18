function [K] = makeOrbitalMatrix(L)
%SUMMARY: This function makes the transformation matrix
% associated with computing the real orbital functions from 
% the complex ones. The dimension of the returned matrix is 
% 2*L + 1. 

I = eye(L);
A = diag((-1).^(1:L));

E = zeros(L,L);
for j = 1:L
E(L - (j - 1),j) = 1;
end

temp = sqrt(1/2)*[1i*I, -1i*A*E;E, A*I];

K = eye(2*L + 1);

ind = 1:(2*L + 1);
ind(L + 1) = [];
K(ind,ind) = temp;

end

