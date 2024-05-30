function [success,r_min,Y,Z,W] = runcap(a,c,q,nu,nubar,Njac,Nalias,Nfft,Ncol,Nrow,Ntail)
% runcap attempts to perform the computer assisted proof
%
% Input:
% a: two-dimensional matrix containing the Fourier coefficients
% (only positive part: cosine series representation)
% c: Wave speed
% q: Frequency
% nu: Weight corresponding to the l^1_nu norm (1x2 array), nu>=1
% nubar: Weight used in the aliassing error estimate
% Njac,Nalias,Nfft,Ncol,Nrow,Ntail: Computional parameters (1x2 arrays)
%
% Output:
% success: true or false depending on whether the attempted proof was successfull
% r_min: Validation radius
% Y,Z,W: Value of the Y,Z,W-bound

disp('Start calculation Y-bound')

Adag = DF_sb_2D(mid(a),mid(c),mid(q),Nfft,Njac);
A_N = intval(inv(mid(Adag)));

% Tools needed for calculating aliasing error
[C,Cprime,Cdoubleprime] = toolsAliasing(a,nubar,Nfft);

% Calculation Y-bound
Y = Ybound(a,A_N,c,q,nu,nubar,Njac,Nalias,Nfft,C); 
X = ['The value of the Y-bound is ', num2str(sup(Y)), '.'];
disp(X)

% Calculation Z-bound
disp('Start calculation Z-bound (note: this might take a while)')
Z = Zbound(a,A_N,c,q,nu,nubar,Njac,Nalias,Nfft,Ncol,Nrow,Ntail,Cprime);
X = ['The value of the Z-bound is ', num2str(sup(Z)), '.'];
disp(X)

% Calculation W-bound
disp('Start calculation W-bound and the proof')
% Calculation W-bound except for multiplication with exp(rstar)
What = Wbound(a,A_N,c,q,nu,nubar,Njac,Nalias,Nfft,Cdoubleprime);

% Rigorous proof (+ computation W bound)
[success,r_min,~,W,r_star] = radiithm(Y,Z,What);
X = ['The value of the W-bound is ',num2str(sup(W)), ' for r_star = ', num2str(sup(r_star)),'.'];
disp(X)

% Display whether the proof was successful or not
if success
    X = ['Proof was successful! The value of r_min is ', num2str(sup(r_min)),'.'];
    disp(X)
else
    disp('Proof was NOT successful')
end

Y=sup(Y);
Z=sup(Z);
W=sup(W);

end