% z = pwr(x,p);
% Equivalent to x.^p, but more efficient for scalar, integer or half-integer p. (x is
% assumed to be non-sparse.)
% Examples (MATLAB R2013a):

% a.^3 is slower than pwr(a,3).
%
% clear, a = .5+rand(1e4); tic, b = a.^3; toc, d = b-a.*a.*a; max(abs(d(:)))
% Elapsed time is 3.693817 seconds.
% ans =
%      8.881784197001252e-16
% clear, a = .5+rand(1e4); tic, b = pwr(a,3); toc, d = b-a.*a.*a; max(abs(d(:)))
% Elapsed time is 0.333767 seconds.
% ans =
%      0

% a.^.5 is slower than pwr(a,.5) (implemented as sqrt(a)).
%
% clear, a = .5+rand(1e4); tic, b = a.^.5; toc, d = b.*b-a; max(abs(d(:)))
% Elapsed time is 0.483987 seconds.
% ans =
%      2.220446049250313e-16
% clear, a = .5+rand(1e4); tic, b = pwr(a,.5); toc, d = b.*b-a; max(abs(d(:)))
% Elapsed time is 0.297585 seconds.
% ans =
%      2.220446049250313e-16

function z = pwr(x,p)
if ~isscalar(p) || 2*p~=real(fix(2*p))
    z = x.^p;
    return
end
% elementwise power x.^p for integer p
if p==0
    z = ones(size(x));
    return
end
if p<0
    x = 1./x;
    p = -p;
end
if p~=fix(p)
    x = sqrt(x);
    p = 2*p;
end
z = x;
b = 2^fix(log2(p));
p = p-b;
while b>1
    z = z.*z;
    b = fix(b/2);
    if p>=b
        z = z.*x;
        p = p-b;
    end
end
