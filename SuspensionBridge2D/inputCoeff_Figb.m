function a_in = inputCoeff_Figb(N)
% inputCoeff_Figb contains the input Fourier coefficients
% corresponding to a numerical approximation with a single peak at the
% bottom
% c = 1.3 and q = (0.05,0.1)
% Note: these are only the Fourier coefficients with positive indices
%
% Input:
% N: Integer specifying how many coefficients should be used. 
% If N(1)>60 or N(2)>20, the Fourier coefficients are padded with zeros 
%
% Output:
% a_in: Input Fourier coefficients for Newton's method (only positive part)

a_in_struct = load('inputCoeff_Figb.mat');
a_in = a_in_struct.a_out;

a_inpad=zeros(N+1); % Pad with zeros
a_inpad(1:size(a_in,1),1:size(a_in,2)) = a_in;
a_in = a_inpad;
a_in = a_in(1:N(1)+1,1:N(2)+1); % For N(1)<60 or N(2)<20
end