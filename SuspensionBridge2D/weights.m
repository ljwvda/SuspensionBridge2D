function [w] = weights(nu,N)
% weights computes the weights corresponding to the l^1_nu norm
%
% Input:
% nu: Weight corresponding to the l^1_nu norm (1x2 array), nu>=1
% N: 1x2 array containing the size-1 of the object you want to compute the
% norm of
%
% Output:
% w: The weights used in l^1_nu, size N(1)+1 by N(2)+1

[k1,k2] = ndgrid(0:N(1),0:N(2));

w = 4*nu(1).^k1.*nu(2).^k2;

% Take the boundaries into account
w(1,:)=2*nu(1).^k1(1,:).*nu(2).^k2(1,:);
w(:,1)=2*nu(1).^k1(:,1).*nu(2).^k2(:,1);
w(1,1) = 1; 

end