function [J,nj] = get_levels(j)

J = -j:j;
J = [ones(numel(J),1)*j, J'];
nj = size(J,1);

end