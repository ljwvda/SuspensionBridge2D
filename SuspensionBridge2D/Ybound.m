function Y = Ybound(a,A_N,c,q,nu,nubar,Njac,Nalias,Nfft,C)
% Ybound calculates the Y-bound necessary for the existence proof
%
% Input:
% a: two-dimensional matrix containing the Fourier coefficients
% (only positive part: cosine series representation)
% A_N: two-dimensional matrix containing the finite part of the 
% approximation of the derivative of F
% c: Wave speed
% q: Frequency
% Njac: 1x2 array holding the number of Jacobian rows and columns 
% nu: Weight corresponding to the l^1_nu norm (1x2 array), nu>=1
% nubar: Weight used in the aliasing error estimate
% We assume that nubar >= nu
% Njac,Nalias,Nfft: 1x2 array holding the number of rows and columns for these N
% C: Constant holding value of C
%
% Output: 
% Y: Value of the Ybound

%% Computation YNalias
% Calculation first summation of YNalias
% We calculate more elements of F than we need for the first summation of
% YNalias. However, the extra elements are used for the second summation
residue_Nalias = F_sb_2D(a,c,q,Nfft,Nalias); 
resmat_Nalias = reshape(residue_Nalias,Nalias+1);
resmat = resmat_Nalias(1:Njac(1)+1,1:Njac(2)+1); 

% To include the enclosure interval we need to calculate epsbar:
epsbar_fgeo = geometric([1/(nubar(1)^(2*Nfft(1))),1/(nubar(2)^(2*Nfft(2)))],[0,0]);
[k1,k2] = ndgrid(0:Nalias(1),0:Nalias(2));
w = nubar(1).^k1.*nubar(2).^k2;
epsbar = w*epsbar_fgeo;

% Include enclosure interval
resmat = resmat + C*epsbar(1:Njac(1)+1,1:Njac(2)+1)*infsup(-1,1); 
residue = resmat(:);
anresmat = reshape(A_N*residue,Njac(1)+1,Njac(2)+1);
first = norml1nu(anresmat,nu);

% Calculation second summation of YNalias
% Computation linear grid to calculate lambda
[k1,k2] = ndgrid(0:Nalias(1),0:Nalias(2));
% Computation lambda
linear = -c^2*k1.^2*q(1)^2+k1.^4*q(1)^4+2*k1.^2.*k2.^2*q(1)^2*q(2)^2+k2.^4*q(2)^4+ones(size(k1));
linearinv = intval('1')./linear;

% Computation weights
omega_short = weights(nu,Nalias);
% We are only dealing with the summation starting at Njac+1 here, 
% therefore we set the rest equal to 0
omega_short(1:Njac(1)+1,1:Njac(2)+1) = intval(zeros(Njac+1)); 

% Summation second term
% Note that we can use resmat_Nalias here as the linear part is zero
% outside I_Njac anyway
second = sum((abs(resmat_Nalias)+C*abs(epsbar)).*omega_short.*linearinv,'all');

% Overall YNalias
YNalias = first+second;

%% Calculation second part of Y
lambdamininv = intval('1')/lambdamin(Nalias,c,q); 
secY = C*lambdamininv*geometric([nu(1)/nubar(1),nu(2)/nubar(2)],Nalias);

%% Combining to obtain overall Y-bound
Y = YNalias + secY;
% Take the upper bound
Y = intval(sup(Y)); 
end