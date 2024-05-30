function [success,r_min,r_max,W,r_star] = radiithm(Y,Z,What)
% radiithm executes the rigorous existence proof
%
% Input:
% Y: Interval value of the Y-bound
% Z: Interval value of the Z-bound
% What: W-bound except for factor exp(r_star)
%
% Output:
% success: Boolean representing whether the proof was successful or not. If
% success=1 the proof was successful, otherwise it was not
% r_min, r_max: Minimum/maximum value of the radius (distance between true
% solution and numerical approximation) for which the proof is successful
% W: Value of W-bound (depending on r_star)
% r_star: A priori estimate radius

% Initializing
success = false;
r_min = 0; r_max = 0;

%% Calculation W-bound
% Taking midpoint of Y- and Z-bound to make an a priori estimate of the
% radius
Ynum = mid(Y);
Znum = mid(Z);
Whatnum = mid(What);

% Radii polynomial
pnum = @(r) 0.5*Whatnum*exp(r)*r^2-(1-Znum)*r+Ynum;

% A priori estimate radius: r_star
options = optimset('TolFun',1e-8,'TolX',1e-8);
r_starnum = fminbnd(@(r) pnum(r),0,1e-1,options); % Find local minimizer
r_star = intval(r_starnum); % Including interval arithmetics

% Calculation W-bound
W = What*exp(r_star);
W = intval(sup(W));

%% Existence proof
discr = (1-Z)^2-2*Y*W; % Discriminant
% Check conditions of the theorem and calculate the minimum and maximum
% radius
if sup(Z)<1 && inf(discr)>0
    r_min = sup((1-Z-sqrt(discr)) / W);
    r_max = inf(min((1-Z) / W, r_star));
    if r_min<r_max
        success = true; % Proof is successful
    end
end

end