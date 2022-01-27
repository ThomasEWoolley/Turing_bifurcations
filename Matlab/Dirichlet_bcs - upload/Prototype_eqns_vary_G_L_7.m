%% Code for cleaning folders and memory
close all;clc;keep pphome;
listing=dir;
Ids=[listing.isdir];
Ids(1:2)=[];

for i=find(Ids)
    rmdir(listing(i+2).name,'s')
end

p=[]; p=stanparam(p); % Setting up basic p structure
p.plot.auxdict={'F', 'G','d','C','max u_1','min u_1'}; % Parameter names (not needed, but nice)
par=           [ 2   ,0,  9,  0.5]; % Parameters used in functions and Jacobian
p.fuha.outfu=@outfn; % Output function to be plotted, names are in above.
p.nc.ilam=2; % Parameter number to continue over
p.plot.bpcmp=5; % Which component of output function to plot.
p.plot.pcmp=1; % Component for plotting
huclean(p); % Set up figures in nice way

p.sw.sfem=-1; % Use OOPDE settings. Otherwise 0/1 for full/preassembled FEM settings
p.nc.neq=2; % Number of equations

p.fuha.sG=@sG; % Functions to be solved
p.sw.jac=0; % Use numeric Jacobian.

p.dim=2; % Spatial dimensions

sw.sym=2; % Symmetry setting for grid 2 = full symmetric (criss-cross), 1=pseudo criss-cross, 0/non less symmetry. Symmetry is need for the bifurcations to be found.
lx=5; ly=lx; nx=10; ny=round(ly*nx/lx); % Space setup and number of elements to be used, space is [-lx,lx]X[-ly,ly].
pde=stanpdeo2D(lx,ly,nx,ny,sw); % PDE setup
p.fuha.cfu=@cf; % Semi-implicit part (in functions and Jacobian)

p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=2*p.np; p.nt=pde.grid.nElements;
p.sw.verb=0; % Verbosity of calculation
p.sol.xi=1/p.nu; % Normalisation weights

p.sw.bifcheck=2; % Calculate bifurcations using eigenvalues, 0 turn off, 1 use LU decomposition
p.sw.spcalc=1; % Calculate eigenvalues (0/1 Eigenvalue computations off/on)
p.nc.neig=min(50,p.np); % Number of eigenvalues to calculate
p.nc.tol=1e-11; % Tolerance for finding branch

p.sol.ds=0.01; % Arclength step
p.nc.dsmin=1e-4; % Arclength minimum steplength
p.nc.dsmax=0.01; % Arclength max steplength
p.nc.dlammax=0.1; % Maximum step in the parameter
p.nc.lammax=3; % Max parameter value
p.nc.lammin=0; % Min parameter value
% p.nc.mu2=0.05; % Threshold for bifcheck=2
p.file.smod=10; % Output every n steps

p.plot.pstyle=1; %0 only plot the FEM mesh; 1 mesh-plot; 2 density-plot
p.sw.para=1; % Continuation variable 1: automatic switching via ? <> p.nc.lamdtol (0: natural parameter.; 2: arclength).


u=ones(p.np,1); v=ones(p.np,1); %ICs
p.u=[u; v; par']; % Complete description

p=oosetfemops(p);
p.nc.nsteps=50;
p=setfn(p,'hom1D',0.01);
p=cont(p,500);
%%


Checked=[];


listing=dir;
L={listing.name};
L={L{[listing.isdir]}};
L1={L{3:end}};
for i=1:length(Checked)
    idx = ~(strcmp([L1{:}], Checked{i}))
    L1=[L1{idx}];
end

while length(L1)>0
    DirName=L1{1};
    
    FileName=dir(['./',DirName,'/bpt*']);
    for i=1:length(FileName)
        try
            p=swibra(DirName,FileName(i).name(1:end-4),[DirName,FileName(i).name(1:end-4)],0.1);
            p.sw.newt=0; % 1: Newton-Cotes
            p.nc.imax=40;
            p.sw.bifcheck=1; % Calculate bifurcations using eigenvalues, 0 turn off, 1 use LU decomposition
            p.sw.foldcheck=1;
            p.sw.spcalc=1; % Calculate eigenvalues (0/1 Eigenvalue computations off/on)
            p.sol.ds=0.001; % Arclength step
            p.nc.tol=1e-8; % Tolerance for residual
            p.nc.dsmin=1e-6; % Arclength minimum steplength
            p.nc.dsmax=0.01;
            p.sw.runpar=0;
            p.plot.pstyle=1;
            p=cont(p,500);
        end
    end
    Checked{1+length(Checked)}=DirName;
    
    listing=dir;
    L={listing.name};
    L={L{[listing.isdir]}};
    L1={L{3:end}};
    for i=1:length(Checked)
        idx = ~(strcmp(L1, Checked{i}));
        L1=sort([L1(idx)]);
    end
    [~,b]=sort(cellfun(@length,L1));
    L1=L1(b);
    
    
end
% %%
% close all
% clc
% plotbra('hom1D')
% plotbra('str2')
% plotbra('str1')
% plotbra('bpt1str1')
% plotbra('bpt2str1')
% plotbra('bpt3str1')
% plotbra('bpt1str2')
% plotbra('bpt2str2')
% plotbra('bpt3str2')
%
% figure
% plotsol('bpt1str1','pt20')



function r=sG(p,u) % pde for chemotaxis model
par=u(p.nu+1:end); %{'a', 'b','d'}
u=u(1:p.nu);n=p.np;u1=u(1:n); u2=u(n+1:2*n); % Seperate parameters and variables
F=par(1);G=par(2);d=par(3);C=par(4);
f1=-u1+u1.^2.*u2-C*(u1-1).^3;
f2=(G+1)+(F-G)*u2-(F+1)*u1.*u2;
f=[f1;f2]; % Nonlinear kinetics

gr=p.pdeo.grid; fem=p.pdeo.fem; % Grid components
K=p.mat.K; N=spd(zeros(n,1));

r=[K N; N d*K]*[u1;u2]-p.mat.M*f; % putting the equation together
end

function A=spd(v)
n=length(v);
A=spdiags(v,0,n,n);
end


function out=outfn(p,u)
% output to bifurcation diagram function
out=[u(p.nu+1:end); % parameters
    max(u(1:p.np));
    min(u(1:p.np))];
end