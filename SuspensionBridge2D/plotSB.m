function plotSB(a_out,q,Nplot,Climit,v)
% plotSB creates a plot of the numerical approximation
% Input:
% a_out: Fourier coefficients that are improved using Newton's method (only positive part)
% q: Frequency
% Nplot: Number of modes used to create the plot
% Climit: Limits of colormap
% v: Angle for viewing the plot

% Fourier transformation and ensuring the solution is symmetric
a_outsym = [a_out flip(a_out(1:end,2:end),2); flip(a_out(2:end,1:end),1) ...
    rot90(a_out(2:end,2:end),2)];
u = ifft2(a_outsym, 'symmetric')*numel(a_outsym);
u2 = u(1:size(a_out,1),1:size(a_out,2));
u2sym = [rot90(u2(2:end,2:end),2) flip(u2(2:end,1:end)); flip(u2(1:end,2:end),2) u2 ];

% Create grid with size of u2sym
y1bold = linspace(-pi/q(1), pi/q(1), size(u2sym,1));
y2bold = linspace(-pi/q(2), pi/q(2), size(u2sym,2));

% Create grid for plotting: use extra modes to make plot more smooth
y1b = linspace(-pi/q(1), pi/q(1), Nplot(1));
y2b = linspace(-pi/q(2), pi/q(2), Nplot(2));
[Xb, Yb] = ndgrid(y1b, y2b);

% Interpolate u2sym onto the new grid
u2sym_interp = interp2(y1bold,y2bold,u2sym',Xb,Yb,'spline'); 

% Create plot
h=surf(Xb,Yb,u2sym_interp);
xlabel('$x_1$','Interpreter','Latex','FontSize',20);
ylabel('$x_2$','Interpreter','Latex','FontSize',20);
zlabel('$\bar{u}$','Interpreter','Latex','FontSize',20);

% Set properties plot
shading interp
lightangle(-25,40)
h.FaceLighting = 'gouraud';
load('cmapEx.mat','cmap');
colormap(cmap);
view(v);
clim(Climit);
xlim('tight');
ylim('tight');

end