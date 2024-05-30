function f = geometric(xi,N)
% geometric calculates f_geo
%
% Input:
% xi: 1x2 array where each element is in [0,1)
% N: 1x2 array
%
% Output:
% f: Value of f_geo 

% Calculate numerator
f_num1 = 2*(xi(2))^(N(2)+1)*(1+xi(1))+2*(xi(1))^(N(1)+1)*(1+xi(2));
f_num2 = -4*(xi(1))^(N(1)+1)*(xi(2))^(N(2)+1);
f_num = f_num1+f_num2;
% Calculate denominator
f_denom = (1-xi(1))*(1-xi(2));

% Calculate f
f = f_num/f_denom;
end