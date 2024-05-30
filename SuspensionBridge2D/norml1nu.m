function [l1nu] = norml1nu(a,nu)
% norml1nu computes the l^1_nu norm of a
%
% Input:
% a: two-dimensional matrix containing the Fourier coefficients (only positive part)
% nu: Weight corresponding to the l^1_nu norm (1x2 array), nu>=1
%
% Output:
% l1nu: The value of the norm

N=size(a)-1;
wmat = weights(nu,N);

l1nu = sum(wmat.*abs(a),'all'); 

end