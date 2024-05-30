function bprimesum = bprimesum(bprime,k,Nout)
% bprimesum sums the bprimes for all combination of indices (+ and -)
%
% Input:
% bprime: Fourier coefficients of the derivative of the nonlinear term
% k: Column we want to compute bprimesum for
% Nout: 1x2 array holding number of rows and columns of output (Njac in practice)
%
% Output:
% bprimesum: Summation of the 4 cases of bprime

% Summation over the 4 different indices (both indices can be positive and negative)
bprimesum = bprime(1+abs((0:Nout(1))-k(1)),1+abs((0:Nout(2))-k(2))) ...
    + bprime(1+abs((0:Nout(1))+k(1)),1+abs((0:Nout(2))-k(2))) ...
    + bprime(1+abs((0:Nout(1))-k(1)),1+abs((0:Nout(2))+k(2))) ...
    + bprime(1+abs((0:Nout(1))+k(1)),1+abs((0:Nout(2))+k(2)));

% Include factor gamma_k/4
if k(1)==0 && k(2)==0
    bprimesum = bprimesum/4;
elseif k(1)==0 ||  k(2)==0
    bprimesum = bprimesum/2;
end

end
