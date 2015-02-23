clear all;
close all;
addpath('./../Assemble/')
addpath('./../FiniteDifferences/')
addpath('./../GeoIntegrals/')
addpath('./../IndexConversion/')
addpath('./../RealComplexConversion/')
addpath('./../VelocityField/')
addpath('./../WignerSymbols/')
tic;
global params
params.nmax=14;
params.nr=20;
params.diffMethod='symmetric2ndOrder';
global velocity
velocity.referenceCase='Gubbins2000';
%% initialization
initializeDiscretization();
initializeVelocity();
%% assembling
CouplingMatrix=assembleCouplingMatrix();
BesselOperator=assembleBesselOperator();
RSquaredMatrix=assembleRSquaredMatrix();
toc