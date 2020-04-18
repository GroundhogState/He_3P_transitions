function [cg] = cg_coef(j1,j2,j,m1,m2)
%SUMMARY: Computes the Clebsch-Gordan matrix element
% for the joint spin-quantum state |j1,m1>|j2,m2>

m = m1 + m2;

cnum = factorial(1 + j1 + j2 + j)*factorial(j1 + m1)*factorial(j1 - m1)*...
       factorial(j + m)*factorial(j - m)*(2*j + 1);

cdenom = factorial(j1 + j2 - j)*factorial(j1 - j2 + j)*...
         factorial(j + j2 - j1)*factorial(j2 + m2)*factorial(j2 - m2);

c = sqrt(cnum/cdenom);

f = @(t) factorial(j1 + j2 - m - t)*factorial(j2 + j - m1 - t)/...
         (factorial(j1 - m1 - t)*factorial(j - m - t)*...
         factorial(j1 + j2 + j + 1 - t));

s = 0;
for k = 0:min(j - m, j1 - m1)
    s = s + (-1)^(j1 - m1 + k)/factorial(k)*f(k);
end

cg = c*s;

end




