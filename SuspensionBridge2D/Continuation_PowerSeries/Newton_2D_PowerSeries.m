function [a_out] = Newton_2D_PowerSeries(x)
% Newton_2D_PowerSeries performs Newton's method to improve Fourier
% coefficients, using the method based on power series. For this method we
% have chosen the cut-off value of the summation of the power series of the
% nonlinearity to be 15, i.e. M=15
%
% Input:
% x: One-dimensional vector containing the Fourier coefficients (a), the
% wave speed (c), the frequency (q), the size of a-1 (n), and the
% cut-off value for the power series of the nonlinearity (M)
%
% Output:
% a_out: Improved Fourier coefficients

tol = 1.e-14; % Tolerance

F = F_sb_2D_PowerSeries(x);
itercount = 0;
while (itercount <= 30) && (max(abs(F)) > tol)
    X = ['Biggest absolute value of F is ',num2str(max(max(abs(F))))];
    display(X)
    J = DF_sb_2D_PowerSeries(x);
    % We do not need the derivative with respect to c here
    J = J(:,1:end-1); 
    a = x(1:end-5);
    a = a(:);
    a = a-J\F; % Newton's method
    x(1:end-5) = a;
    F = F_sb_2D_PowerSeries(x);
    itercount = itercount+1;
end
if max(abs(F)) > tol          %%% The corrector did not converge %%%
    disp('Failure! Initial guess is not good enough for the Newton algorithm to converge')
    a_out = a_in;
    return
end

a_out = x;
end