
function gauss=custom_gaussian_function_2d(fun_size, sigma, center, theta, offset, amp)
% gaussian_function_2d - Compute a two dimensional gaussian with widths,rotation & offsets
% based on code by Thomas Dideriksen %https://au.mathworks.com/matlabcentral/fileexchange/9556-custom-2d-gauss
%   -improved vectorization
%   -center is excatly at center of matrix even if matrix is even sized
% Matlabs inbilt gaussian function fspecial('gaussian', [30 30], 4) has a the same sd in each dirn
% which makes it inapropriate to use as a blur kernel for 'squished' 2d plots.
% does not inculde rotation
% the middle of the output matrix coresponds to the (0,0) of the function
%
% Syntax:  [data,import_opts]=import_data(import_opts)
%
% Inputs:
%    filter_size	- vector[1x2], of the output matrix
%    sigma          - vector[1x2], width of the gaussian in units of the filter_size
%
% Outputs:
%    gauss      - matrix, sampled gaussian
%
%
% Example: 
%   imagesc(Gaussian_filter2d([100 100],[10 10]))
%   set(gca,'dataAspectRatio',[1 1 1])

% Other m-files required: test_custom_gaussian_function_2d
% Also See: none
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%   - add optional input 'normalization',['peak','sum','none']
%
% Author: Bryce Henson
% email: Bryce.m.henson+gaussian_function_2d@gmail.com
% Last revision:2019-03-07

%------------- BEGIN CODE --------------
[xvals,yvals] = meshgrid((1:fun_size(1))-((fun_size(1)-1)/2)-1,...
    (1:fun_size(2))-((fun_size(2)-1)/2)-1);
xc      = center(1);
yc      = center(2);
theta   = (theta/180)*pi;
xm      = (xvals-xc).*cos(theta) - (yvals-yc).*sin(theta);
ym      = (xvals-xc).*sin(theta) + (yvals-yc).*cos(theta);
gauss   = exp(-((xm./sigma(1)).^2 + (ym./sigma(2)).^2)./2);

if amp~=0 || ~isnan(amp) %replace this with optional argument
    gauss= gauss/sum(gauss(:));
    gauss=gauss*amp;
end
gauss=gauss+offset;

end
