% Jeremy & Wes

winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight



dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);

% Variables
tSim = 600e-15; % simulation running time
f = 230e12; % frequency of plane wave
lambda = c_c/f; % wavelength of plane wave

% Region setup
xMax{1} = 20e-6;    % max x coordinate in region
nx{1} = 200;    % number of nodes in region in x
ny{1} = 0.75*nx{1}; % number of nodes in region in y

% Characteristics of region
Reg.n = 1; % number of different regions

mu{1} = ones(nx{1},ny{1})*c_mu_0; % permeability of region

epi{1} = ones(nx{1},ny{1})*c_eps_0; % permittivity of region 1
epi{1}(125:150,55:95)= c_eps_0*11.3; % permittivity of region 2 within region 1

sigma{1} = zeros(nx{1},ny{1});  % conductivity of region 1
sigmaH{1} = zeros(nx{1},ny{1}); % magnetic conductivity of region 1

% Region measurements
dx = xMax{1}/nx{1}; % x step length
dt = 0.25*dx/c_c; % time step length
nSteps = round(tSim/dt*2); % number of steps in simulation
yMax = ny{1}*dx; % maximum y coordinate in region
nsteps_lamda = lambda/dx; % number of wavelength steps thru region

% Plotting settings
movie = 1; % play movie yes or no
Plot.off = 0; % dont plot yes or no
Plot.pl = 0;
Plot.ori = '13'; % subplot orientation
Plot.N = 100; % plot N number of steps
Plot.MaxEz = 1.1; % maximum Ez value
Plot.MaxH = Plot.MaxEz/c_eta_0; % maximum H value
Plot.pv = [0 0 90]; % custom plot viewing angle
Plot.reglim = [0 xMax{1} 0 yMax]; % region limits

% Boundary Conditions
bc{1}.NumS = 2; % number of sources
bc{1}.s(1).xpos = nx{1}/(4) + 1; % soft source starting position
bc{1}.s(2).xpos = 3*nx{1}/(4) - 1; % soft source starting position
bc{1}.s(1).type = 'ss'; % type of soft source
bc{1}.s(2).type = 'ss'; % type of soft source
bc{1}.s(1).fct = @PlaneWaveBC; % use PlaneWaveBC function to determine E and H values
bc{1}.s(2).fct = @PlaneWaveBC; % use PlaneWaveBC function to determine E and H values
% mag = -1/c_eta_0;
mag = 1; % magnitude of H field
phi = 0; % phase shift
omega = f*2*pi; % wave velocity
betap = 0; % propagation constant
t0 = 30e-15; % time shift
% st = 15e-15; % 
st = -0.05;
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'p'}; % soft source parameters s - spherical
bc{1}.s(2).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'}; % soft source parameters p - plane

Plot.y0 = round(y0/dx);

bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

pml{1}.width = 20 * spatialFactor;
pml{2}.width = 0.5 * pml{1}.width;
pml{1}.m = 3.5;
pml{2}.m = 0.5 * pml{1}.m;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg

% 2ci) pml is the inclusion, but by commenting it out, we remove variables used
% in RunYeeReg.m

% 2cii) bc structure defines all boundary conditions for the system

% 2ciii) bc{1}.s(1) gives the parameters for soft source #1

% 2civ) bc{1}.xm/xp/ym/yp affects ###

% 3b) st does something with the left side of the region

% 3c) increasing the frequency too high causes the plane wave to not
% penetrate through the grating, decreasing the frequency gives more
% prominent plane wave peaks

% 4b) added second source, swapped both sources between spherical and plane
% waves

