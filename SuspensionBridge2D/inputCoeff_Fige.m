function a_in = inputCoeff_Fige(N)
% inputCoeff_smallc contains the input Fourier coefficients
% corresponding to a numerical approximation with two peaks at the bottom
% c = 0.9 and q = (0.1,0.1)
% Note: these are only the Fourier coefficients with positive indices
%
% Input:
% N: Integer specifying how many coefficients should be used. 
% If N(1)>60 or N(2)>60, the Fourier coefficients are padded with zeros 
%
% Output:
% a_in: Input Fourier coefficients for Newton's method

a_in_struct = load('inputCoeff_Fige.mat');
a_in = a_in_struct.a_out;

a_inpad=zeros(N+1); % Pad with zeros
a_inpad(1:size(a_in,1),1:size(a_in,2)) = a_in;
a_in = a_inpad;
a_in = a_in(1:N(1)+1,1:N(2)+1); % For N(1)<60 or N(2)<60
end