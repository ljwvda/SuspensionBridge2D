function [matrixnorm] = matrixnorml1nu(A,nu,n)
% norm of matrix as operator from l^1_nu to l^1_nu
%
% Input:
% A: two-dimensional matrix containing the Fourier coefficients (only positive part)
% nu: Weight corresponding to the l^1_nu norm (1x2 array), nu>=1
% n: size of Fourier coefficients -1
%
% Output:
% matrixnorm: The value of the norm

wmat = weights(nu,n); % Calculate weights
v = 1./wmat;

% Calculate operator norm
matrixnorm = max((wmat(:)'*abs(A)).*v(:)'); 

end