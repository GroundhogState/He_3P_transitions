function test_custom_gaussian_function_2d
mat_size=[1e3 1e3];
sigma=[10 10];
tic
filt1=custom_gaussian_function_2d(mat_size, sigma, [0,0], 0, 0, 1);
time_custom=toc;
tic
filt2=fspecial('gaussian', mat_size, sigma(1)) ;
time_inbuilt=toc;
filt2=filt2/sum(filt2(:));
diff=filt1-filt2;
max_diff=max(diff(:));
imagesc(diff)

amp_tol=1e-9;
fprintf('testing custom against inbuilt matlab 2d gaussian\n')
logic_str = {'FAIL', 'pass'};
fprintf('INFO: custom runtime       : %g\n',time_custom)
fprintf('INFO: inbuilt runtime      : %g\n',time_inbuilt)
fprintf('INFO: max difference       : %g\n',max_diff)

fprintf('TEST: difference within tolerance: %s\n',logic_str{1+(max_diff<amp_tol)})

end
