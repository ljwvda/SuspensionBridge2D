function [J] = DF_sb_2D(a,c,q,Nfft,Nout)
% DF_SB_2D calculates the derivative of F
%
% Input:
% a: Two-dimensional matrix containing the Fourier coefficients (only positive part)
% c: Wave speed
% q: Frequency
% Nfft:  1x2 array holding the number of rows and columns for this N
% Nout: 1x2 array holding number of rows and columns of output J
%
% Output: 
% J: Derivative as a 2D object

[k1,k2] = ndgrid(0:Nout(1),0:Nout(2)); % Grid containing the indices

Nxy=(Nout(1)+1)*(Nout(2)+1);

% Linear part
linear = -c^2*k1.^2*q(1)^2+k1.^4*q(1)^4+2*k1.^2.*k2.^2*q(1)^2*q(2)^2+k2.^4*q(2)^4+ones(size(k1));
Jlinear = diag(linear(:));

% Nonlinear part
Jnonlinear = nonlinearityMatrix(a,Nfft,Nout);
Jnonlinear = reshape(Jnonlinear,Nxy,Nxy);

% Combining linear and nonlinear part
J = Jlinear+Jnonlinear;
end