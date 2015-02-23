function [A,B]=getEigProb(Remag,nmax,nr)
%%
% getEigProb -- returns sparse matrices A and B to solve the eigenvalue problem
%
%     A*v = \lambda B*v
%
% whereas B is a diagonal matrices with non zero entries on the diagonal. Hence,
% the matrix B is easy to invert.
%
% function call:
% 
%     [A,B]=getEigProb(Remag,nmax,nr);
%
% If the input arguments are not specified, default values are used. That is
%
%   Remag=100, nmax=14, nr=25.
%
%%
if not(nargin==3)
	Remag=100;
	nmax=14;
	nr=25;
elseif nargin==3
else
	error('Please call this funciton with zero or three inputs.')
end
if exist('OCTAVE_VERSION')
	strVersion=strjoin({'GNU Octave',version});
else
	strVersion=strjoin({'MATLAB',version});
end
%------------------------------------------------------------------------------%
%% start of the programm -- welcome message
%------------------------------------------------------------------------------%
tic;
fprintf(['------------------------------------------------------------\n',...
'------------------------------------------------------------\n',...
'       ____        ____               __ \n',...
'      / __ )__  __/ / /___ __________/ / \n',...
'     / __  / / / / / / __ `/ ___/ __  /  \n',...
'    / /_/ / /_/ / / / /_/ / /  / /_/ /   \n',...
'   /_____/\\__,_/_/_/\\__,_/_/   \\__,_/    \n',...
'     ______     ____                    \n',...
'    / ____/__  / / /___ ___  ____ _____ \n',...
'   / / __/ _ \\/ / / __ `__ \\/ __ `/ __ \\ \n',...
'  / /_/ /  __/ / / / / / / / /_/ / / / / \n',...
'  \\____/\\___/_/_/_/ /_/ /_/\\__,_/_/ /_/  \n',...
'          ______          __   \n',...
'         / ____/___  ____/ /__  \n',...
'        / /   / __ \\/ __  / _ \\ \n',...
'       / /___/ /_/ / /_/ /  __/ \n',...
'       \\____/\\____/\\__,_/\\___/  \n\n',...
'       authors --- Sebastian Glane\n',...
'       version --- 0.9\n',...
' executable in --- GNU Octave 3.8.1\n',...
'               --- MATLAB R2013a (8.1.0.604)\n',...
'   your system --- %s\n',...
'------------------------------------------------------------\n'],strVersion);
displayOutput();
%------------------------------------------------------------------------------%
%% specifying input parameters
%------------------------------------------------------------------------------%
%% global parameter structure for discretization
global params
params.nmax=nmax;
params.nr=nr;
params.diffMethod='symmetric2ndOrder';
%% global parameter structure for velocity field
global velocity
% velocity.referenceCase='ToroidalTestCase';
velocity.referenceCase='Gubbins2000';
velocity.M=-0.08;
velocity.D=0.5;
velocity.p=3;
%------------------------------------------------------------------------------%
%% initialization phase
%------------------------------------------------------------------------------%
fprintf(['------------------------------------------------------------\n',...
'------- initialization phase -------------------------------\n',...
'------------------------------------------------------------\n'])
displayOutput();
ticid=tic;
%% adding path for discretization
addpath('./lib/FiniteDifferences/')
initializeDiscretization();
%% adding pathes for velocity initialization
addpath('./lib/RealComplexConversion/')
addpath('./lib/SpecialFunctions/')
addpath('./lib/VelocityField/')
initializeVelocity();
fprintf(['------------------------------------------------------------\n',...
'------- end of initialization phase --- timing: %8.5f ---\n',...
'------------------------------------------------------------\n'],toc(ticid));
displayOutput();
%------------------------------------------------------------------------------%
%% assembling phase
%------------------------------------------------------------------------------%
fprintf(['------------------------------------------------------------\n',...
'------- assembly phase -------------------------------------\n',...
'------------------------------------------------------------\n']);
displayOutput();
ticid=tic;
%% adding pathes for assemble
addpath('./lib/Assemble/')
addpath('./lib/IndexConversion/')
addpath('./lib/GeoIntegrals/')
addpath('./lib/WignerSymbols/')
CouplingMatrix=assembleCouplingMatrix();
BesselOperator=assembleBesselOperator();
RSquaredMatrix=assembleRSquaredMatrix();
fprintf(['------------------------------------------------------------\n',...
'------- end of assembly phase --- timing: %8.5f ---------\n',...
'------------------------------------------------------------\n'],toc(ticid));
displayOutput();
%------------------------------------------------------------------------------%
%% function output
%------------------------------------------------------------------------------%
A=BesselOperator+Remag*CouplingMatrix;
B=RSquaredMatrix;
end