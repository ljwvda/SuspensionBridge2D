function [C,Cprime,Cdoubleprime] = toolsAliasing(a,nubar,Nfft)
% toolsAliasing calculates constants that appear in the aliasing error
% see Section 3.1 and appendix B of the paper
%
% Input:
% a: Fourier coefficients in 2D (only positive part)
% nubar: Weight used in the aliasing error estimate
% We assume that nubar >= nu
% Nfft: 1x2 array holding the number of rows and columns for this N
%
% Output:
% C,Cprime,Cdoubleprime: Constants C,C' and C'' (note that we do not have 
% the hat as we are working with cosine series)

% rhobar needed for calculation C, C' and C''
rhobar(1) = log(nubar(1));
rhobar(2) = log(nubar(2));

% delta are the elements of a small rectangle with sides
% [0,pi/Nfft(1)], [0,pi/Nfft(2)]
delta(1) = infsup(0,sup(intval('pi')/(Nfft(1)))); 
delta(2) = infsup(0,sup(intval('pi')/(Nfft(2))));
% Grid of elements in J_Nfft in the paper
[n1a,n2a] = ndgrid(-Nfft(1):Nfft(1)-1,-Nfft(2):Nfft(2)-1);
% Power of exponential in abar_n^{rhobar,delta} in paper
power = n1a*rhobar(1)+n2a*rhobar(2)+1i*(n1a*delta(1)+n2a*delta(2)); 

% Pad Fourier coefficients and make them symmetric
apad = intval(zeros(Nfft+1)); 
apad(1:size(a,1),1:size(a,2))=a; 
a_sympad = [rot90(apad(2:end,2:end),2) flip(apad(2:end,1:end-1),1); ...
    flip(apad(1:end-1,2:end),2) apad(1:end-1,1:end-1)]; 
% Calculation abar_n^{rhobar,delta} in paper
abar_delta = a_sympad.*exp(power); 
% Note that the size of abar_delta is a power of 2, 
% which is necessary for altfftn
% FFT to go back to ubar
ubar_delta = altifftn(abar_delta)*(size(abar_delta,1)*size(abar_delta,2)); 

% Computation C, C' and C''
% Use G
C = 1/(4*Nfft(1)*Nfft(2))*sum(abs(exp(ubar_delta)-ubar_delta-1), 'all');
C = intval(sup(C));

% Use G'
Cprime = 1/(4*Nfft(1)*Nfft(2))*sum(abs(exp(ubar_delta)-1), 'all');
Cprime = intval(sup(Cprime));

% Use G''
Cdoubleprime = 1/(4*Nfft(1)*Nfft(2))*sum(abs(exp(ubar_delta)), 'all');
Cdoubleprime = intval(sup(Cdoubleprime));

end