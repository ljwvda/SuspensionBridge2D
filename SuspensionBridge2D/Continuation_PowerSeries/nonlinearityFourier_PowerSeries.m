function [vp,vpfull] = nonlinearityFourier_PowerSeries(v,M,first)
% nonlinearityFourier_PowerSeries applies Fourier transformations and
% corresponding shifts and calculates nonlinear part
% Also works for higher dimensions
%
% Input:
% v: Fourier coefficients cosine series
% M: Cut-off value for the power series of the nonlinearity
% First: Optional input. Boolean indicating whether the first term of the
% summation has to be included
%
% Output:
% vp: has the same size as v , coefficients of nonlinear part
% vpfull: Contains all the positive modes. If M=0, vpfull has size(v)

if ~exist('first','var')
    % If third parameter does not exist, default it to something
    first = 0;
end

sz = size(v);
dim=length(sz);

N = (sz-1)/2;

if M==0
    % special case of power=zero
    vp=0*v;
    for k=1:dim
        s{k}=N(k)+1;
    end
    vp(s{1:dim})=1;
    vpfull=vp;
    return
end

Mp = 2.^nextpow2(2*M*N+1);

v1 = altzeros(Mp(1:dim),v(1));
M2=Mp/2;
M2(Mp==1)=0; % Take care of dimensions with just one element (N=0)
s1 = M2 + 1 - N;
s2 = M2 + 1 + N;
s1f = M2 + 1 - M*N;
s2f = M2 + 1 + M*N;
for j=1:dim
    s{j}=s1(j):s2(j);
    sf{j}=s1f(j):s2f(j);
    fl{j}=[M2(j)+1:Mp(j),1:M2(j)];
end

v1(s{1:dim})=v;
v2 = v1(fl{1:dim}); %fftshift
w = altfftn(v2);

% Determine if first term should be added
if first == 1
    wp = w;
    for i = 2:(M-1)
        wp = wp + w.^i/factorial(i);
    end
else
    
    wp = 0;
    for i = 2:M
        wp = wp + w.^i/factorial(i);
    end
end

vp2 = altifftn(wp);
vp1 = vp2(fl{1:dim}); %fftshift

vp = vp1(s{1:dim});
vpfull = vp1(sf{1:dim});

return



