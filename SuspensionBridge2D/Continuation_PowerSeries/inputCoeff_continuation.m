function a_in = inputCoeff_continuation(nx,ny)
% inputCoeff_continuation contains the input Fourier coefficients
% corresponding to a numerical approximation of the Suspension Bridge
% equation for c=1.4, q=(0.1,0.05)
% Note: these are only the Fourier coefficients with positive indices
%
% Input:
% nx: Integer specifying how many rows should be used. 
% ny: Integer specifying how many columns should be used. 
% If nx>43 or ny>43, the Fourier coefficients are padded with zeros 
%
% Output:
% a_in: Input Fourier coefficients for Newton's method

a_inpad=zeros(nx,ny); % Pad with zeros
a_in_struct = load('inputCoeff_continuation.mat');
a_in = a_in_struct.a_out;
a_inpad(1:size(a_in,1),1:size(a_in,2)) = a_in;
a_in = a_inpad;
a_in = a_in(1:nx,1:ny); % For nx<43 or ny<43
end