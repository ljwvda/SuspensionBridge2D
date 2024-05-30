%%% Main script used for obtaining the example CAP (Theorem 1.1)/Figure a in the paper %%%
% Script can easily be generalized to create other CAPs, just change
% the input (starting at line 17)
% Note: INTerval LABoratory is necessary to obtain the (rigorous) proofs

tic % Measure computation time
close all; clear;

% Add intlab
if ~exist('intval', 'file')
    addpath('~/Intlab_V12'); % Adjust this to your computer
    startintlab;
end

%% Initialization

% Initialize different N's
Ngal = [130 130];
Njac = [65 65];
Nalias = [400 400];
Nfft = [1024 1024];
Ncol = [300 300];
Nrow = [800 800];
Ntail = [140 140];

% Input parameters to obtain example solution
a_in = inputCoeff_Example(Ngal);

c = intval('1.1'); % Include interval arithmetic
q(1) = intval('pi')/intval('50');
q(2) = intval('pi')/intval('40');

nu(1) = intval('1.0000001'); % nu >= 1
nu(2) = intval('1.0000001'); 

a_out = Newton_2D(a_in,mid(c),mid(q),Nfft); % Apply Newton's method
a = intval(a_out); % Include interval arithmetic

% Value of nubar (should be > nu)
nubar(1) = intval('1.1');
nubar(2) = intval('1.1');
%% CAP
[success,r_min,Y,Z,W] = runcap(a,c,q,nu,nubar,Njac,Nalias,Nfft,Ncol,Nrow,Ntail);

%% Plotting
Nplot = [299 299]; % Number of modes used in plotting
Climit = [-10 1.23393]; % Limits of colormap
viewAngle = [-1.21,13.75]; % Angle for viewing the plot
% viewAngle = [-10.5 13]; % Angle for creating Figure a
plotSB(a_out,mid(q),Nplot,Climit,viewAngle)
xlim([-50 50]); % Ensure that the whole domain is visible on the axes
ylim([-40 40]); 

toc