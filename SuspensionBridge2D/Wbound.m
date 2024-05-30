function What = Wbound(a,A_N,c,q,nu,nubar,Njac,Nalias,Nfft,Cdoubleprime)
% Wbound calculates W-bound up to a factor exp(rstar)
%
% Input:
% a: two-dimensional matrix containing the Fourier coefficients
% (only positive part: cosine series representation)
% A_N: two-dimensional matrix containing the finite part of the 
% approximation of the derivative of F
% c: Wave speed
% q: Frequency
% nu: Weight corresponding to the l^1_nu norm (1x2 array), nu>=1
% nubar: Weight used in the aliasing error estimate
% We assume that nubar >= nu
% Njac, Nalias, Nfft: 1x2 array holding the number of rows and columns for these N
% Cdoubleprime: Constant holding value of C''
%
% Output:
% What (almost the W bound): 
% product of l^1_nu norm of matrix A and estimate for the l^1_nu norm of exp(a)

%% Calculate estimate for norm of exponential of a
% Pad a with zeros so it becomes the size we need
apad = intval(zeros(Nalias+1)); 
apad(1:size(a,1),1:size(a,2))=a; 

% Computation bbarFFT_doubleprime
bbarFFT_doubleprime = nonlinearity(apad,2,Nfft,Nalias); 

% Weights
omega = weights(nu,Nalias);

% Calculation epsbar
epsbar_fgeo = geometric([1/(nubar(1)^(2*Nfft(1))),1/(nubar(2)^(2*Nfft(2)))],[0,0]);
[k1small,k2small] = ndgrid(0:Nalias(1),0:Nalias(2));
w = nubar(1).^k1small.*nubar(2).^k2small;
epsbar = w*epsbar_fgeo;

% Summation
first = sum((abs(bbarFFT_doubleprime)+Cdoubleprime*abs(epsbar)).*omega,'all');

geom = geometric([nu(1)/nubar(1) nu(2)/nubar(2)],Nalias); % Calculation f_geo   
% Estimate norm exponential of a
normexpa = first+Cdoubleprime*geom; 

%% Calculate norm of A
% Computation norm of A, divided in finite (n<=Njac) and infinite part
firstA = matrixnorml1nu(A_N,nu,Njac); 
% For the infinite part we only need the smallest eigenvalue, which yields
% the largest inverse eigenvalue
secondA = 1/lambdamin(Njac,c,q);
Anorm = max(firstA,secondA); 

%% Combine
What = Anorm*normexpa;
% Take the upper bound
What = intval(sup(What)); 
end