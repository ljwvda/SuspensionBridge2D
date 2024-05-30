function curve=cont(f,df,x0,q,n,M,nsteps,step0,stepmin,stepmax,tol,silent)
%CONT finds the solution curve to a system of equations
%
%  X=CONT(F,DF,X0,NSTEPS,STEP0,STEPMIN,STEPMAX,TOL,SILENT)
%  X=CONT(F,DF,X0,NSTEPS,STEP0)
%  computes the solution curve X
%  for N-1 equations F=0 in N variables,
%  starting from an initial solution X0,
%  using pseudo-arclength continuation.
%
%  Here X0 should be a column vector of length N and
%  the function F should send a vector of length N
%  to a column vector of length N-1,
%  i.e. F : R^N -> R^(N-1).
%  DF should be the (N-1) x N derivative matrix of F.
%  All vectors are assumed to be column vectors.
%
%  The output X is a matrix with per column a solution of F=0.
%
%  F should be a function_handle. It can be either a handle to
%  an anonymous or to a function file (but not an inline function).
%  If DF=[] then the derivative will be computed numerically.
%
%  In total a maximum of NSTEPS steps are taken,
%  starting with stepsize STEP0 and terminating
%  if the stepsize becomes smaller than STEPMIN.
%  The maximum stepsize is STEPMAX.
%  The initial guess X0 should be no further than STEPMAX from a real zero.
%  The sign of STEP0 determines the direction of the first step
%
%  Default values for STEPMIN and STEPMAX are 10^-8 and 1.
%  Stepsize are measured in terms of L^2 norms, so for large problems
%  the user may want to scale the max/min step size appropriately.
%
%  TOL determines the tolerance in determining zeros (default 10^-14).
%  By default a plot is produced of the first two
%  components of the solution curve,
%  except if the dimension equals 3, when a curve in 3D is plotted.
%
%  If SILENT=1 then no plot is produced and no warning messages are displayed.
%

% default values
if ~exist('stepmax') || isempty(stepmax)
    stepmax=1;
end
if ~exist('stepmin') || isempty(stepmin)
    stepmin=1e-8;
end
if ~exist('tol') || isempty(tol)
    tol=1e-14;
end
if ~exist('silent')
    silent=0;
end

if isempty(df)
    % numerical derivative
    df=@(x) numjac(@(t,y) f([y;q(:);n(:);M]),0,x,f(x),eps);
end

% initialization
dim=length(x0);
y0=reshape(x0,dim,1);
curve=[];
niter=100;
step=sign(step0)*max(min(abs(step0),stepmax),stepmin);
k=0;
%relative amount predictor is allowed to be off
correctiemax=1/10;

% compute first point
gr=df([y0;q(:);n(:);M]);
nullgr=null(gr)';
[ynew,dummy,iters]=multinewton(@(x) [f([x;q(:);n(:);M]);nullgr*(x-y0)],...
    @(x) [df([x;q(:);n(:);M]);nullgr],y0,tol,niter,1);

% check first point
if iters > niter || norm(ynew-y0)>stepmax
    if silent~=1
        fprintf(['could not find a zero near the given starting point.\n'])
    end
    return
end

% start loop
while k<=nsteps & abs(step)>=stepmin
    %update
    y=ynew;
    k=k+1
    curve=[curve,y];
    gr=df([y;q(:);n(:);M]);
    nullgrold=nullgr;
    nullgr=null(gr)';
    
    % check that the gradient did not flip
    if nullgrold*nullgr'<0
        nullgr=-nullgr;
    end
    result=0;
    
    % iterate until successful
    while abs(step)>=stepmin & result~=1
        % predictor
        z=y+step*nullgr';
        % corrector
        [ynew,dummy,iters]=multinewton(@(x) [f([x;q(:);n(:);M]);nullgr*(x-z)],...
            @(x) [df([x;q(:);n(:);M]);nullgr],z,tol,niter,1);
        % check newtonoutput is OK
        if norm(z-ynew)<correctiemax*abs(step) & iters<niter
            result=1;
            % iteration succesfull: adapt stepsize
            % increase stepsize if iters<5
            % decrease stepsize if iters>5
            step=sign(step)*max(stepmin,min(stepmax,abs(step)*2^((4-iters)/3)));
        else
            % corrector too large (or nonconvergence):
            % decrease stepsize and try again
            step=step*2^(-1/3);
        end
    end
end

% check untimely end
if k<nsteps && silent~=1
    fprintf(['Continuation terminated after %i steps ', ...
        'due to lack of convergence.\n'],k-1);
end

% plot
if silent~=1 && k>1
    if dim==3
        plot3(curve(1,:)',curve(2,:)',curve(3,:)','.-')
    else
    % Plot sqrt(c) against l^1_nu norm for nu=1
        plot(sqrt(curve(end,:))', sum(abs(curve(1:end-1,:)))','.-') %added this line
        hold on
        xlabel('$c$','Interpreter','Latex','FontSize',20)
        ylabel('$\|\bar{a}\|_\nu$','Interpreter','Latex','FontSize',20)
    end
end

return


