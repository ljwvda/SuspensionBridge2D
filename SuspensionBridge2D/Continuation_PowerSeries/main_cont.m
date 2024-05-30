%%% This is the main script you need to run %%%

% This code uses pseudo-arclength continuation to continue the solution of
% the suspension bridge equation in the parameter c (wave speed)

close all; clear;

%% Initialization
c2_in=1.4^2;
q(1)=0.1; q(2)=0.05;

a_in_start = inputCoeff_continuation(43,43);

n(1) = size(a_in_start,1)-1;
n(2) = size(a_in_start,2)-1;

M = 15; % For Newton's method wehave chosen the cut-off value of the 
% summation of the power series of the nonlinearity to be 15
x = [a_in_start(:);c2_in;q(:);n(:);M];

%% Newton's method
% Improve numerical approximation
x_out = Newton_2D_PowerSeries(x);

%% Continuation (plotting)
% Make sure the input is in the correct manner
x_out = x_out(:);
a = x_out(1:end-6);

tol=1e-12; % Tolerance
% Continuation for M=14, plot is made inside cont.m
M = 14; % Cut-off value for the power series of the nonlinearity
cont(@F_sb_2D_PowerSeries,@DF_sb_2D_PowerSeries,[a;c2_in],q,n,M,50,-0.001,[],[],tol);
hold on
% Continuation for M=15, plot is made inside cont.m
M = 15; % Cut-off value for the power series of the nonlinearity
cont(@F_sb_2D_PowerSeries,@DF_sb_2D_PowerSeries,[a;c2_in],q,n,M,50,-0.001,[],[],tol);
legend('M=14','M=15','Interpreter','Latex','FontSize',20) % Add legend
