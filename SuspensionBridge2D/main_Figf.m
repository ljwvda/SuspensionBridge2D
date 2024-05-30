%%% Main script used for obtaining the CAP in Figure 3f for large c in the paper %%%
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
Ngal = [80 40];
Njac = [80 60];
Nalias = [200 150];
Nfft = [512 512];
Ncol = [220 170];
Nrow = [1000 1000];
Ntail = [100 100];

% Input parameters to obtain example solution
a_in = inputCoeff_Figf(Ngal);

c = intval('1.4'); % Include interval arithmetic
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
Nplot = [299 299]; % Number of modes used in plotting
Climit = [-10 1.23393]; % Limits of colormap
viewAngle = [-10.5 13]; % Angle for viewing the plot
plotSB(a_out,mid(q),Nplot,Climit,viewAngle)

toc