function Z = Zbound(a,A_N,c,q,nu,nubar,Njac,Nalias,Nfft,Ncol,Nrow,Ntail,Cprime)
% Zbound calculates the Z-bound necessary for the existence proof
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
% Njac,Nalias,Nfft,Ncol,Nrow,Ntail: 1x2 array holding the number of rows and columns for these N
% where we assume Njac<=Ntail<=Ncol and Njac<=Nrow
% Cprime: Constant holding value of C'
%
% Output:
% Z: Value of the Z-bound

%% Computing variables needed for computing Z per column
% We calculate these outside the for-loops to speed up the computation
% Calculation bbar_prime
% Computation for n in Nalias

% Pad a with zeros so it becomes the size we need
apad = intval(zeros(Nalias+1));
apad(1:size(a,1),1:size(a,2))=a;

% Calculation bbarFFT_prime
bbarFFT_prime_small = nonlinearity(apad,1,Nfft,Nalias);

% Calculation epsbar
epsbar_fgeo = geometric([1/(nubar(1)^(2*Nfft(1))),1/(nubar(2)^(2*Nfft(2)))],[0,0]);
[k1small,k2small] = ndgrid(0:Nalias(1),0:Nalias(2));
w = nubar(1).^k1small.*nubar(2).^k2small;
epsbar = w*epsbar_fgeo;

% Calculation bbar' for small n
bbar_prime_small = bbarFFT_prime_small + Cprime*epsbar*infsup(-1,1);

% Estimate bbar' for large n
[k1,k2] = ndgrid(0:Ncol(1)+Nrow(1),0:Ncol(2)+Nrow(2));
bbar_prime = infsup(-1,1)*Cprime./(nubar(1).^abs(k1).*nubar(2).^abs(k2));

% Combine
bbar_prime(1:Nalias(1)+1,1:Nalias(2)+1) = bbar_prime_small;

% Calculation rest of the tools
% Linear grid
[k1,k2] = ndgrid(0:Njac(1),0:Njac(2));
linear = -c^2*k1.^2*q(1)^2+k1.^4*q(1)^4+2*k1.^2.*k2.^2*q(1)^2*q(2)^2+k2.^4*q(2)^4+ones(size(k1));

% Weights
omega = weights(nu,max([Nalias;Ncol;Nrow]));

% Reshape A_N to 4D object
AN4D = reshape(A_N,Njac(1)+1,Njac(2)+1,Njac(1)+1,Njac(2)+1);

% Calculation linear part of F for size Nrow
[k1large,k2large] = ndgrid(0:Nrow(1),0:Nrow(2));
linearlarge = -c^2*k1large.^2*q(1)^2+k1large.^4*...
    q(1)^4+2*k1large.^2.*k2large.^2*q(1)^2*q(2)^2+k2large.^4*q(2)^4+ones(size(k1large));
% Inverting and setting the part that was already included for small n to zero
linearinvlarge = intval('1')./linearlarge;
linearinvlarge(1:Njac(1)+1,1:Njac(2)+1) = intval(zeros(Njac+1));

geomNrow = geometric([nu(1)/nubar(1) nu(2)/nubar(2)],Nrow); % Calculation f_geo
% Smallest value of lambda outside I^+_Nrow
lambdamin_Nrow = lambdamin(Nrow,c,q);

%% For "small" columns, we compute the Zbound per column and store the
% maximum value
maxZcol = intval(0); % Initialize maximum value
for c1 = 0:Ncol(1)
    for c2 = 0:Ncol(2)
        [Zcol,Zcolall] = Zbound_col(A_N,Njac,Nrow,nu,nubar,[c1 c2],linear,...
            omega,bbar_prime,Cprime,AN4D,linearinvlarge,geomNrow,lambdamin_Nrow);
        if sup(Zcol)>sup(maxZcol)
            kcol(1)=c1;
            kcol(2)=c2;
            Zcolmax=Zcolall;
        end
        maxZcol = intval(max(Zcol,maxZcol)); % Store maximum value
        if c1==ceil(Ncol(1)/2) && c2==ceil(Ncol(2)/2) % Check progress
            disp('Halfway through the for-loops for computing the Z-bound per column')
        end
    end
end

%% For "large" columns, we estimate the Zbound
%% Estimate for small n
muNtail = muhat(Ntail,Ncol,nu,nubar);
muNjac = muNtail(1:Njac(1)+1,1:Njac(2)+1);
ANmu = abs(A_N)*muNjac(:);
ANmu = reshape(ANmu,Njac+1);
Zest1 = Cprime*sum(ANmu.*omega((1:Njac(1)+1),(1:Njac(2)+1)),'all');

%% Estimate for intermediate n
Zest2 = Cprime*sum(muNtail.*linearinvlarge((1:Ntail(1)+1),(1:Ntail(2)+1)).*omega((1:Ntail(1)+1),(1:Ntail(2)+1)),'all');

%% Estimate for large n
% Part of bprime necessary for the estimate
bbar_prime = bbar_prime(1:Nalias(1)+1,1:Nalias(2)+1);

est_largen_1 = sum(abs(bbar_prime).*omega((1:Nalias(1)+1),(1:Nalias(2)+1)),'all');
geomNalias = geometric([nu(1)/nubar(1) nu(2)/nubar(2)],Nalias); % Calculation f_geo
est_largen_2 = Cprime*geomNalias;

lambdaminNtail= lambdamin(Ntail,c,q);
Zest3 = intval('1')/lambdaminNtail*(est_largen_1+est_largen_2);

%% Combining
Zest = Zest1 + Zest2 + Zest3;

%% Overall Z-bound
Z = max(maxZcol,Zest);
Z = intval(sup(Z));

%% output for debugging and optimization purposes
disp(['The maximum in Zcol is attained at kcol = [ ',num2str(kcol(1)),' , ',num2str(kcol(2)),' ]']);
disp(['The corresponding value of Zcol1 is ',num2str(sup(Zcolmax(1)))]);
disp(['The corresponding value of Zcol2 is ',num2str(sup(Zcolmax(2)))]);
disp(['The corresponding value of Zcol3 is ',num2str(sup(Zcolmax(3)))]);
disp(['The uniform bound Zest1 = ',num2str(sup(Zest1))]);
disp(['The uniform bound Zest2 = ',num2str(sup(Zest2))]);
disp(['The uniform bound Zest3 = ',num2str(sup(Zest3))]);
disp(['The value of Zcol is ',num2str(sup(maxZcol))]);
disp(['The value of Zest is ',num2str(sup(Zest))]);
end
