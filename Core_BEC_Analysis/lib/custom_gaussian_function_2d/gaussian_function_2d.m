
function gauss=legacy_gaussian_function_2d(filter_size, sigma)
%LEGACY WRAPER to keep functionality of gaussian_function_2d intact
gauss=custom_gaussian_function_2d(filter_size, sigma, [0,0], 0, 0, 1);
end
