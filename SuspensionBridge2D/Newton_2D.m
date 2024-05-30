function [a_out] = Newton_2D(a_in,c,q,Nfft)
% Newton_2D performs Newton's method to improve input
%
% Input:
% a_in: Two-dimensional matrix containing the Fourier coefficients (only positive part)
% c: Wave speed
% q: Frequency
% Nfft: 1x2 array holding the number of rows and columns for this N
%
% Output:
% a_out: Improved Fourier coefficients

N = size(a_in)-1;
tol = 1.e-14; % tolerance for Newton's method
x = a_in;  
F = F_sb_2D(x,c,q,Nfft,N);
iter = 0; % Iteration counter

% Repeatedly applying Newton's method
disp('Starting Newton iterations');
while (iter <= 30) && (max(abs(F)) > tol)
    X = ['Biggest absolute value of F is ',num2str(max(abs(F)))];
    disp(X); % Show biggest absolute value of F to check performance	
    J = DF_sb_2D(x,c,q,Nfft,N);
    x = x(:);
    x = x-J\F; 
    x = reshape(x,N+1);
    F = F_sb_2D(x,c,q,Nfft,N);
    iter = iter+1;                            
end

% No convergence
if max(abs(F)) > tol          
  disp('Failure! Initial guess is not good enough for the Newton algorithm to converge')  
  a_out = a_in;
  return
else
  disp('Newton iterations have converged');
end  

a_out = x;

end