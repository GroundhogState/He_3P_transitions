function [J,C] = clebsch_gordan(j1,j2)
%SUMMARY: Computes the possible joint states given the spin angular
% moment quantum numbers j1, j2 in the basis of Clebsch-Gordan 
% coefficients. The function returns the resulting states in the two-column 
% matrix 'J' as given by the relation | j1 - j2 | <= j <= j1 + j2. 
% The matrix 'C' contains the Glebsch-Gordan coefficients. 

jmin = abs(j1 - j2);
jmax = j1 + j2;

increment = jmin:jmax;

J = []; Nj = [];

for k = 1:numel(increment)
    
    [jtemp,n] = get_levels(increment(k));
    Nj = [Nj; n];
    J = [J;jtemp];
    
end

[J1,n1] = get_levels(j1);
[J2,n2] = get_levels(j2);

N = sum(Nj);
C = zeros(N,n1*n2);

for m = 1:N
    for m1 = 1:n1
        for m2 = 1:n2
            if (J1(m1,2) + J2(m2,2)) == J(m,2)
               C(m,m1 + n1*(m2 - 1)) = cg_coef(J1(m1,1),J2(m2,1),J(m,1),J1(m1,2),J2(m2,2)); 
            end
        end
    end
end

end

