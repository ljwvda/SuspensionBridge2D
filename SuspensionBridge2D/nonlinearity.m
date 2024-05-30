function [vshort,vF] = nonlinearity(v,deriv,Nfft,Nout)
% nonlinearity calculates the nonlinear part of F in terms of Fourier
% coefficients for cosine series. This code ensures that we can work with 
% cosine series instead of Fourier series. The true Fourier transformation
% is done in nonlinearityFourier.m
% Code also works for higher dimensions
%
% Input:
% v: Fourier coefficients cosine series
% deriv: 0 means we consider g(u)=exp(u)-u-1, 1 means that we consider the
% derivative of the nonlinear part, which is g'(u)=exp(u)-1 and 2 means
% that we consider the second derivative, i.e. g''(u)=exp(u)
% Nfft: 1x2 array holding the number of rows and columns for this N. It
% determines the size of the Fourier transform of the nonlinear part. The
% Fourier transformations are done in nonlinearityFourier.m, which outputs 
% only the part we actually need
% Nout: 1x2 array determining the number of rows and columns for the
% outputs (see Output)
%
% Output:
% vshort: has size Nout+1
% vF: contains all modes, not just the positive ones. Its size is
% 2*size(vshort)-1, so 2*Nout+1

N = size(v)-1;
dim=length(N);

for j=1:dim
  s{j}=1:N(j)+1;
  sflip{j}=abs(-N(j):N(j))+1;
  sfull{j}=Nout(j)+1:2*Nout(j)+1;
end

vv=v(sflip{1:dim});
% applies Fourier transformations and corresponding shifts and calculates
% nonlinear part
vF=nonlinearityFourier(vv,deriv,Nfft,Nout); 
vF=real(vF);

vshort = vF(sfull{1:dim});

end
