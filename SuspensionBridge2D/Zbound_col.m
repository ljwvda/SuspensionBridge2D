function [Zcol,Zcolall] = Zbound_col(A_N,Njac,Nrow,nu,nubar,col,linear,omega,...
    bbar_prime,Cprime,AN4D,linearinvlarge,geom,lambdamin_Nrow)
% Zbound_col calculates Z-bound for each column (up to Ncol)
%
% Input:
% A_N: two-dimensional matrix containing the finite part of the
% approximation of the derivative of F, size is (n1 x n2) x (n1 x n2)
% Njac, Nrow: 1x2 array holding the number of rows and columns for these N
% nu: Weight corresponding to the l^1_nu norm (1x2 array), nu>=1
% nubar: Weight used in the aliasing error estimate
% We assume that nubar >= nu
% col: Indicates for which 2D-column we evaluate the Zbound, this should be
% strictly positive
% linear,omega,bbar_prime,Cprime,AN4D,mid_bprime,linearinvlarge,geom,lambdamin_Nrow:
% tools calculated in toolsZ.m
%
% Output:
% Zcol: Zbound for column col
% Zcolall: The three components [Zcol1,Zcol2,Zcol3] that sum up to Zcol

%% summation n<=N

% For col<=Njac (first term), we need standard unit vectors
if col(1)<=Njac(1) && col(2)<=Njac(2)
    ID1 = intval(eye(Njac(1)+1)); % For creating unit vector
    evec1 = ID1(:,col(1)+1);
    ID2 = intval(eye(Njac(2)+1));
    evec2 = ID2(:,col(2)+1);
end

% The first term is only non-zero for small columns
if col(1)<=Njac(1) && col(2)<=Njac(2)
    kronecker = evec1*evec2';
    term1 = kronecker - AN4D(:,:,col(1)+1,col(2)+1)*linear(col(1)+1,col(2)+1);
else
    term1 = intval(zeros(Njac+1));
end

% Calculation second term (n<=Njac)
% Summation over the 4 different indices (both indices can be positive and negative)
bprime4 = bprimesum(bbar_prime,col,Njac);
term2 = A_N*bprime4(:);
term2 = reshape(term2,Njac+1);

smalln_num = sum(abs(term1-term2).*omega((1:Njac(1)+1),(1:Njac(2)+1)),'all');
Zcol1 = smalln_num/omega(col(1)+1,col(2)+1);

%% large n: n>N
% First term
bprime4_large = abs(bprimesum(bbar_prime,col,Nrow)); 
largen_first_num = sum(linearinvlarge.*bprime4_large.*omega((1:Nrow(1)+1),(1:Nrow(2)+1)),'all');
Zcol2 = largen_first_num/omega(col(1)+1,col(2)+1);

% Estimate
fract_num = Cprime*nubar(1)^col(1)*nubar(2)^col(2);
fract_denom = lambdamin_Nrow*nu(1)^col(1)*nu(2)^col(2);
Zcol3 = fract_num/fract_denom*geom;

%% Combining everything

Zcol=Zcol1+Zcol2+Zcol3;
Zcolall=[Zcol1,Zcol2,Zcol3];

end