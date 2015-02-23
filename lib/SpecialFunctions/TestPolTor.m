clear all
close all;
clc;
addpath('./../vtkExport/')
dimensions=[3,50,50];
radius=1.0;
Data.dimensions=dimensions;
Data.radius=radius;
[R,Theta,Phi]=getSphericalGrid(radius,dimensions);
Ind=[3 5];
%------------------------------------------------------------------------------%
[Field]=getToroidalFieldReal(Ind,[R Theta Phi]);
FieldCos=Field(:,:,1);
FieldSin=Field(:,:,2);
%------------------------------------------------------------------------------%
[CartesianFieldCos]=SphericalField2Cartesian(FieldCos,[R Theta Phi]);
[CartesianFieldSin]=SphericalField2Cartesian(FieldSin,[R Theta Phi]);
%------------------------------------------------------------------------------%
Data.Vectors=cell(2,2);
Data.Vectors{1,1}=CartesianFieldCos;
Data.Vectors{1,2}='ToroidalFieldCos';
Data.Vectors{2,1}=CartesianFieldSin;
Data.Vectors{2,2}='ToroidalFieldSin';
%------------------------------------------------------------------------------%
% fpath='./ToroidalDump.vtk';
% ftitle='Toroidal Vector Field'
% writeData2vtk(Data,fpath,ftitle)
% 
%------------------------------------------------------------------------------%
% [Field]=getPoloidalFieldReal(Ind,[R Theta Phi]);
% FieldCos=Field(:,:,1);
% FieldSin=Field(:,:,2);
% 
%------------------------------------------------------------------------------%
% [CartesianFieldCos]=SphericalField2Cartesian(FieldCos,[R Theta Phi]);
% [CartesianFieldSin]=SphericalField2Cartesian(FieldSin,[R Theta Phi]);
% 
%------------------------------------------------------------------------------%
% Data.Vectors=cell(2,2);
% Data.Vectors{1,1}=CartesianFieldCos;
% Data.Vectors{1,2}='PoloidalFieldCos';
% Data.Vectors{2,1}=CartesianFieldSin;
% Data.Vectors{2,2}='PoloidalFieldSin';
% 
%------------------------------------------------------------------------------%
% fpath='./PoloidalDump.vtk';
% ftitle='Poloidal Vector Field'
% writeData2vtk(Data,fpath,ftitle)
%------------------------------------------------------------------------------%