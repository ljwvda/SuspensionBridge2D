function [J] = DF_sb_2D_PowerSeries(x)
% DF_sb_2D_M14 computes DF corresponding to the partial derivatives of the 
% zero-finding problem of the suspension bridge equation using the method 
% based on power series
% Taking M=14 terms in the series into account
% Note that we now also compute the derivative with respect to c^2
%
% Input:
% x: One-dimensional vector containing the Fourier coefficients (a), the
% wave speed squared (c^2), the frequency (q), the size of a-1 (n), and the
% cut-off value for the power series of the nonlinearity (M)
%
% Output:
% J: Two-dimensional matrix containing the derivatives. The last column is
% the derivative with respect to c

% Split x into the different variables
a=x(1:end-6);
c2=x(end-5);
q(1) = x(end-4);q(2) = x(end-3);
n(1)=x(end-2);n(2)=x(end-1);
M=x(end);

Nxy=(n(1)+1)*(n(2)+1);

[k1,k2] = ndgrid(0:n(1),0:n(2));

a=reshape(a,n+1);

% Linear term
linear = -c2*k1.^2*q(1)^2+k1.^4*q(1)^4+2*k1.^2.*k2.^2*q(1)^2*q(2)^2+k2.^4*q(2)^4+1; 
Jlinear = diag(linear(:));

% Nonlinear term
nonlinear = nonlinearityMatrix_PowerSeries(a,n,M); 

% Reshape as a matrix
Jnonlinear=reshape(nonlinear,Nxy,Nxy);

% Combining linear and nonlinear part
DFa = Jlinear + Jnonlinear;

% Adding the derivative with respect to c^2
DFc2 = -k1.^2*q(1)^2.*a;
DFc2=DFc2(:);

J = [DFa DFc2];

end