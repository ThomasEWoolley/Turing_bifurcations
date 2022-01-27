%% Code for cleaning folders and memory
close all;clc;keep pphome;
listing=dir;
Ids=[listing.isdir];
Ids(1:2)=[];

for i=find(Ids)
    rmdir(listing(i+2).name,'s')
end

p=[]; p=stanparam(p); % Setting up basic p structure
p.plot.auxdict={'L','max u_1','min u_1'}; % Parameter names (not needed, but nice)
par=           [ 0.04]; % Parameters used in functions and Jacobian
p.fuha.outfu=@outfn; % Output function to be plotted, names are in above.
p.nc.ilam=1; % Parameter number to continue over
p.plot.bpcmp=2; % Which component of output function to plot.
p.plot.pcmp=2; % Component for plotting
huclean(p); % Set up figures in nice way

p.sw.sfem=-1; % Use OOPDE settings. Otherwise 0/1 for full/preassembled FEM settings
p.nc.neq=2; % Number of equations

p.fuha.sG=@sG; % Functions to be solved
p.sw.jac=0; % Use numeric Jacobian.

p.dim=1; % Spatial dimensions

lx=1; nx=50;% Space setup and number of elements to be used, space is [-lx,lx]X[-ly,ly].
pde=stanpdeo1D(0.5,0.004); % PDE setup
% bc=pde.grid.robinBC('1','1'); 
p.fuha.cfu=@cf; % Semi-implicit part (in functions and Jacobian)

p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=2*p.np; p.nt=pde.grid.nElements;
p.sw.verb=0; % Verbosity of calculation
p.sol.xi=1/p.nu; % Normalisation weights

p.sw.bifcheck=2; % Calculate bifurcations using eigenvalues, 0 turn off, 1 use LU decomposition
p.sw.spcalc=1; % Calculate eigenvalues (0/1 Eigenvalue computations off/on)
p.nc.neig=min(50,p.np); % Number of eigenvalues to calculate
p.nc.tol=1e-11; % Tolerance for finding branch

p.sol.ds=0.001; % Arclength step
p.nc.dsmin=1e-7; % Arclength minimum steplength
p.nc.dsmax=0.01; % Arclength max steplength
p.nc.dlammax=0.1; % Maximum step in the parameter
p.nc.lammax=0.5; % Max parameter value
p.nc.lammin=0; % Min parameter value
% p.nc.mu2=0.05; % Threshold for bifcheck=2
p.file.smod=10; % Output every n steps

p.plot.pstyle=1; %0 only plot the FEM mesh; 1 mesh-plot; 2 density-plot
p.sw.para=1; % Continuation variable 1: automatic switching via ? <> p.nc.lamdtol (0: natural parameter.; 2: arclength).


u=2*ones(p.np,1); v=0.75*ones(p.np,1); %ICs
p.u=[u; v; par']; % Complete description

p=oosetfemops(p);
p.nc.nsteps=50;
p=setfn(p,'hom1D',0.01);
p=cont(p,100);
%%
 huclean(p); % Set up figures in nice way
p=swibra('hom1D','bpt1','1D1',-0.001);
p.plot.bpcmp=2;
p.plot.pcmp=1; % Component for plotting
p.nc.dsmax=0.1; % Arclength max steplength
p.nc.lammax=2; % Max parameter value
p=cont(p,2000);
%%
huclean(p); % Set up figures in nice way
p=swibra('hom1D','bpt1','1D2',0.001);
p.plot.pcmp=1; % Component for plotting
p.plot.bpcmp=3;
p.nc.dsmax=0.01; % Arclength max steplength
p.nc.lammax=2; % Max parameter value
p=cont(p,2000);

%%
% huclean(p); % Set up figures in nice way
p=swibra('hom1D','bpt1','1D2',-0.001);
p.plot.pcmp=1; % Component for plotting
p.nc.dsmax=0.1; % Arclength max steplength
p.nc.lammax=2; % Max parameter value
p=cont(p,200);

%% check near bif point
% huclean(p); % Set up figures in nice way
p=swibra('hom1D','bpt1','Near bif',0.001);
p.plot.pcmp=1; % Component for plotting
p.nc.dsmax=0.01; % Arclength max steplength
p.nc.lammax=2; % Max parameter value
p=cont(p,20);

%%
% p=swibra('1D1','bpt1','bp1_from1D1',-0.001);
p0=qswibra('1D1','bpt1');
p=seltau(p0,1,'bp1_from1D1',1);
p.plot.pcmp=1; % Component for plotting
p.nc.dsmax=0.1; % Arclength max steplength
p.nc.lammax=3; % Max parameter value
p=cont(p,200);



function r=sG(p,u) % pde for Schnak
par=u(p.nu+1:end); %{'a', 'b','d'}
u=u(1:p.nu);n=p.np;u1=u(1:n); u2=u(n+1:2*n); % Seperate parameters and variables
L=par(1);
f1=1-2*u1+u1.^2.*u2;
f2=3-u1.^2.*u2;
f=[f1;f2]; % Nonlinear kinetics

gr=p.pdeo.grid; fem=p.pdeo.fem; % Grid components
K=p.mat.K; N=spd(zeros(n,1));

Qu=p.mat.Qu; Qv=p.mat.Qv; bcu=p.mat.Gu; bcv=p.mat.Gv;
r=[1e-3*K/L^2 N; N 1e-1*K/L^2]*[u1;u2]-p.mat.M*f; % putting the equation together
end

function A=spd(v)
n=length(v);
A=spdiags(v,0,n,n);
end


function out=outfn(p,u)
% output to bifurcation diagram function
% u=u(1:p.nu);n=p.np;u1=u(1:n); u2=u(n+1:2*n); % Seperate parameters and variables
out=[u(p.nu+1:end); % parameters
    max(u(1:p.np));
    min(u(1:p.np))];
end