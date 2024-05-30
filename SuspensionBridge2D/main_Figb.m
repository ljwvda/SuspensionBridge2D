%%% Main script used for obtaining the CAP in Figure 3b in the paper %%%
% Note: INTerval LABoratory is necessary to obtain the proofs

tic % Measure computation time
close all; clear;

% Add intlab
if ~exist('intval', 'file')
    addpath('~/Intlab_V12'); % Adjust this to your computer
    startintlab;
end

%% Initialization

% Initialize different N's
Ngal = [60 20];
Njac = [60 20];
Nalias = [150 150];
Nfft = [1024 1024];
Ncol = [165 160];
Nrow = [380 380];
Ntail  = [80,30];

% Input parameters to obtain example solution
a_in = inputCoeff_Figb(Ngal);

c = intval('1.3'); % Include interval arithmetic
q(1) = intval('0.05');
q(2) = intval('0.1');

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
Nplot = [399 399]; % Number of modes used in plotting
Climit = [-10 1.23393]; % Limits of colormap
viewAngle = [-10.5 13]; % Angle for viewing the plot
plotSB(a_out,mid(q),Nplot,Climit,viewAngle)

toc