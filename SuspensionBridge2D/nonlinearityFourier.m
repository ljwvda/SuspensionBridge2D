function [vF] = nonlinearityFourier(v,deriv,Nfft,Nout)
% nonlinearityFourier applies Fourier transformations and corresponding shifts
% and calculates nonlinear part
% Code also works for higher dimensions
%
% Input:
% v: Fourier coefficients
% deriv: 0 means we consider g(u)=exp(u)-u-1, 1 means that we consider the
% derivative of the nonlinear part, which is g'(u)=exp(u)-1 and 2 means
% that we consider the second derivative, i.e., g''(u)=exp(u)
% Nfft: 1x2 array holding the number of rows and columns for this N. It
% determines how much zero padding is used before we apply Fourier transformations
% and it determines the size of the Fourier transformation
% Nout: 1x2 array determining the number of rows and columns for the
% output (see Output)
%
% Output:
% vF: Fourier coefficients of nonlinear part, which has size 2*Nout+1

sz = size(v);
dim=length(sz);

N = (sz-1)/2;

Nfft2 = 2 * Nfft; 

v1 = altzeros(Nfft2(1:dim),v(1)); 

s1 = Nfft + 1 - N;
s2 = Nfft + 1 + N;

for j=1:dim
    s{j}=s1(j):s2(j);
    fl{j}=[Nfft(j)+1:Nfft2(j),1:Nfft(j)];
end

v1(s{1:dim})=v;
v2 = v1(fl{1:dim}); %fftshift
w = altifftn(v2)*prod(size(v2));

if deriv == 0 % G
    wp=exp(w)-w-1;
elseif deriv == 1 % G'
    wp=exp(w)-1;
elseif deriv == 2 % G''
    wp=exp(w);
end

vF2 = altfftn(wp)/prod(size(wp));
vF1 = vF2(fl{1:dim}); %fftshift

% vF1 has size Nfft2
% For vF we only take out a part of size 2*Nout+1, therefore we create sout
s1out = Nfft + 1 - Nout;
s2out = Nfft + 1 + Nout;
for j=1:dim
    sout{j}=s1out(j):s2out(j);
end

vF = vF1(sout{1:dim}); 

return



