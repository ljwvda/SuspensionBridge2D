function [F] = F_sb_2D(a,c,q,Nfft,Nout)
% F_SB_2D calculates F
%
% Input:
% a: Two-dimensional matrix containing the Fourier coefficients (only positive part)
% c: Wave speed
% q: Frequency
% Nfft: 1x2 array holding the number of rows and columns for this N
% Nout: 1x2 array holding number of rows and columns of output F
%
% Output: 
% F: F as a column vector

[k1,k2] = ndgrid(0:Nout(1),0:Nout(2)); % Grid containing the indices

% Linear part
linear = -c^2*k1.^2*q(1)^2+k1.^4*q(1)^4+2*k1.^2.*k2.^2*q(1)^2*q(2)^2+k2.^4*q(2)^4+ones(size(k1));

if isintval(a) % Use intervals if input is interval
    apad = intval(zeros(Nout+1));
else
    apad = zeros(Nout+1);
end
apad(1:size(a,1),1:size(a,2)) = a;
apad = apad(1:Nout(1)+1,1:Nout(2)+1); 

% Nonlinear part
nonlinear = nonlinearity(a,0,Nfft,Nout);

% Combining linear and nonlinear part
F = linear.*apad+nonlinear;

F=F(:); % Reshape F
end