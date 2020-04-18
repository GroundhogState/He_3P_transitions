function [Yr,Yc,xyz] = spherical_harmonics(radius, npts, L)
%SUMMARY: This function takes as arguments 'radius', the radius of the sphere where 
% the spherical harmonic functions takes their values. 'npts' is the sqare
% root of the total number of points that is used to generate the spherical
% point cloud. The point cloud is the unit sphere scaled by 'radius'. 
% 'L' is the orbital angular momentum associated with the spherical
% harmonics. For a given 'L', the function returns the complex and real spherical 
% harmonics associated with 'L' and the magnetic quantum number 'M'
% associated with 'L', where |M| <= L runs over the signed integers. 
% The function returns the cell structures 'Yr' and 'Yc', the real and
% complex wave functions, respectively. The spherical harmonics are ordered
% in the cell structure from beginning to end as
% -M, -M + 1, -M + 2, .. 0, 1, 2, ... M such that 
% there will be totally 2*L + 1 entries in each cell structure. 

%External dependencies: 'makeOrbitalMatrix.m', 'delta_l.m'
%Internal dependencies: 'linspace.m', 'kron.m', 'sin.m','cos.m', 'cell.m', 
%                       'sqrt.m', 'factorial.m', 'zeros.m', 'size.m', ++ ...                        
% The list of internal dependencies is probably not complete. 



N = npts;
theta = linspace(0,pi,N)';
phi = linspace(0,2*pi,N)';

thetphi = [kron(ones(N,1),theta),kron(phi,ones(N,1))];
xyz = radius*[sin(thetphi(:,1)).*cos(thetphi(:,2)), sin(thetphi(:,1)).*sin(thetphi(:,2)), cos(thetphi(:,1))];

principal_n = 2*L + 1;
Yc = cell(1,principal_n);

counter = 0;
for M = -L:L

counter = counter + 1;    
dl = delta_l(L,M);

Alm = sqrt((2*L + 1)*factorial(L - M)*factorial(L + M))/sqrt(4*pi);
Ylm = zeros(size(xyz,1),1);

for j = 1:size(dl,1)

c = factorial(dl(j,1))*factorial(dl(j,2))*factorial(dl(j,3));
pqs = dl(j,:);
Ylm = Ylm + Alm/c * ( -( xyz(:,1) + 1i*xyz(:,2) )/2 ).^(pqs(1)) ...
           .* ( (xyz(:,1) - 1i*xyz(:,2))/2 ).^(pqs(2)).*xyz(:,3).^(pqs(3));

end
Ylm = Ylm/(radius^L);

Yc{1,counter} = Ylm;
end

K = makeOrbitalMatrix(L);

Yr = cell(1,principal_n);
for i = 1:principal_n
    temp = zeros(size(Yc{1,i}));
    for j = 1:principal_n
        temp = temp + K(i,j)*Yc{1,j};
    end
    Yr{1,i} = real(temp);
end
       
end

