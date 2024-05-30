function lambda_min = lambdamin(N,c,q)
% lambdamin calculates smallest value of lambda outside I^+_N
%
% Input: 
% N: 1x2 array holding the number of rows and columns for this N
% c: Wave speed
% q: Frequency
%
% Output:
% lambda_min: Smallest value of lambda outside I^+_N

if N(1)<= min(c/q(1)) % Check requirement
    error('N(1) is chosen too small. Value of lambda_min might be incorrect.')
end

% Define lambda's as grid
[k1,k2] = ndgrid(0:N(1)+1,0:N(2)+1);
linear = -c^2*k1.^2*q(1)^2+k1.^4*q(1)^4+2*k1.^2.*k2.^2*q(1)^2*q(2)^2+k2.^4*q(2)^4+ones(size(k1));

% First value is calculated outside the loop as initialization
lambda_aux = linear(1,N(2)+2); 

% Calculate smallest value of lambda outside I^+_N
for i=1:ceil(sup(c/q(1)))
    if linear(i+1,N(2)+2)<lambda_aux
        lambda_aux = linear(i+1,N(2)+2);
    end
end

lambda_min = min(lambda_aux, linear(N(1)+2,1));

end