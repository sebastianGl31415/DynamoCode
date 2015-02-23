clear all;
close all;
clc;
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
params.nmax=14;
params.nr=50;
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
% eigenvalue computation
%------------------------------------------------------------------------------%
fprintf(['------------------------------------------------------------\n',...
'------- eigenvalue computation -----------------------------\n',...
'------------------------------------------------------------\n']);
return
neig=5;
Remag=[50:7.5:150];
fprintf('       number of eigenvalues to be found --- %d\n',neig);
displayOutput();
Lambda=zeros(neig,length(Remag));
for k=1:length(Remag)
	[V,D,fl]=eigs(BesselOperator+Remag(k)*CouplingMatrix,...
		RSquaredMatrix,neig,'lr');
	Lambda(:,k)=diag(D);
	if fl==0
		fprintf('       Remag --- %4.0f / convergence\n',Remag(k));
	else
		fprintf('       Remag --- %4.0f / no convergence\n',Remag(k));
	end
	displayOutput();
end
fprintf(['------------------------------------------------------------\n',...
'------- end of eigenvalue computation ----------------------\n',...
'------------------------------------------------------------\n']);
%------------------------------------------------------------------------------%
% post processing
%------------------------------------------------------------------------------%
fprintf(['------------------------------------------------------------\n',...
'------- post processing ------------------------------------\n',...
'------------------------------------------------------------\n']);
strM=sprintf('%3.2f',velocity.M);
delpos=strfind(strM,'0.');
strM=strcat(strM(1:delpos-1),strM(delpos+2:end));
strD=sprintf('%3.2f',velocity.D);
delpos=strfind(strD,'0.');
strD=strcat(strD(1:delpos-1),strD(delpos+2:end));
fname=strcat(velocity.referenceCase,'_Nmax',num2str(params.nmax),...
				'_Nr',num2str(params.nr),'_M',strM,'_D',strD,'.mat');
if exist('OCTAVE_VERSION')
		save('-mat-binary',fname,'neig','Lambda','Remag');
else
		save(fname,'neig','Lambda','Remag');
end
fprintf(' written data to  --- %s\n',fname);
fprintf(['------------------------------------------------------------\n',...
'------- end of post processing -----------------------------\n',...
'------------------------------------------------------------\n']);