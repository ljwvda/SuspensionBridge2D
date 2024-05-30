function [vpshort,vpfull,vpextended] = nonlinearity_PowerSeries(v,M,first)
% nonlinearity_PowerSeries calculates the nonlinear part of F in terms of
% Fourier coefficients for cosine series
% Also works for higher dimensions
%
% Input:
% v: Fourier coefficients cosine series
% M: Cut-off value for the power series of the nonlinearity
% First: Optional input. Boolean indicating whether the first term of the
% summation has to be included
%
% Output:
% vpshort: has the same size as v 
% vpfull: contains all the positive modes. If M=0, vpfull has size(v)
% vpextended: contains all modes, not just the positive ones

 if ~exist('first','var')
     % If third parameter does not exist, default it to something
      first = 0;
 end

N = size(v)-1;
dim=length(N);

p0=max(1,M); % Takes care of p=0 case

for j=1:dim
  s{j}=1:N(j)+1;
  sflip{j}=abs(-N(j):N(j))+1;
  sfull{j}=p0*N(j)+1:2*p0*N(j)+1;
end

vv=v(sflip{1:dim});
[~,vpextended]=nonlinearityFourier_PowerSeries(vv,M,first);
vpextended=real(vpextended);

vpfull = vpextended(sfull{1:dim});
vpshort = vpfull(s{1:dim});

return
