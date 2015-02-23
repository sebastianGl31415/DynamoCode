function []=writeData2vtk(Data,fpath,filetitle)
% writeField2Vtk -- writes the data contained in given structure Data to 
%     vtk-file fpath with the title filetitle. The data must be given in 
%     spherical coordinates with number of points in each coordinate 
%     [nr,ntheta,nphi] specified in Data.dimensions and radius given in 
%     Data.radius. The data to write must be stores in the structure Data in the
%     following way:
%       -- scalars: Data.Scalars{:,1:2} is cell array containing vectors of 
%                   length prod(Data.dimensions) and a string for the name.
%       -- vectors: Data.Vectors{:,1:2} is cell array containing matrices of 
%                  size ( prod(Data.dimensions) , 3 ) and a string for the name.
%       -- tensors: Data.Tensors{:,1:2} is cell array containing fields of size
%                   ( prod(Data.dimensions) , 3 , 3 ) and a string for the name.
%
%  []=writeData2vtk(Data,fpath,filetitle)
%
% author: Sebastian Glane

	%% input check
	if not(ischar(fpath))
	error('2nd input is not a string.');
	elseif not(ischar(filetitle))
		error('3rd input is not a string.');
	elseif not(isstruct(Data))
		error('1st input is not a structure.');
	end
	%% getting sphere radius
	if not(isfield(Data,'radius'))
		error('1st input does not contains field radius.')
	else
		radius=Data.radius;
	end
	%% getting point numbers
	if not(isfield(Data,'dimensions'))
		error('1st input does not contains field dimensions.')
	else
		dimensions=Data.dimensions;
	end
	%% getting geometry
	[R,Theta,Phi]=getSphericalGrid(radius,dimensions);
	[Points,delRows]=getCartesianPoints(R,Theta,Phi);
	nPoints=length(Points);
	clear R Theta Phi
	%% removing doubling points from data
	Data=cleanDataset(Data,delRows);
	%% getting Cells and CellTypes
	[Cells]=getCells(dimensions);
	[Cells,CellTypes]=getCellTypes(Cells);
	%% number of celltypes and cells
	nCellTypes=length(Cells);
	nCells=length(CellTypes);
	%% counting number of elements in Cells
	cellListSize=0;
	for k=1:nCellTypes
		cellListSize=cellListSize+numel(Cells{k});
	end
	%% user info
	fprintf('Writing vtk-file: %s ',fpath)
	%% open file
	fid = fopen(fpath,'w');
	if fid==-1
		error('Cannot open the file to write.');
	end
	%% writing vtk-header
	fprintf(fid,'# vtk DataFile Version 2.0 \n');
	fprintf(fid,'%s\n',filetitle);
	fprintf(fid,'BINARY\n');
	fprintf(fid,'DATASET UNSTRUCTURED_GRID \n');
	fprintf(fid,'POINTS %d float\n',nPoints);
	%% writing points
	fwrite(fid,Points','float','b');
	%% writing cells
	fprintf(fid,'\nCELLS %d %d\n',nCells,cellListSize);
	for k=1:nCellTypes
		fwrite(fid,int32(Cells{k}'),'int32','b');
	end
	fprintf(fid,'\nCELL_TYPES %d\n',nCells);
	fwrite(fid,cast(CellTypes,'int32'),'int32','b');
	%% writing data
	fprintf(fid,'\nPOINT_DATA %d\n',nPoints);
	if isfield(Data,'Scalars')
		for k=1:size(Data.Scalars,1)
			if not(isempty(Data.Scalars{k,1}))
				writeData=transpose(Data.Scalars{k,1});
				fprintf(fid,'SCALARS %s float\n',Data.Scalars{k,2});
				fprintf(fid,'LOOKUP_TABLE default %d\n',length(Data.Scalars{k,1}));
				fwrite(fid,writeData,'float','b');
			end
		end
	end
	if isfield(Data,'Vectors')
		for k=1:size(Data.Vectors,1)
			if not(isempty(Data.Vectors{k,1}))
				writeData=transpose(Data.Vectors{k,1});
				fprintf(fid,'\nVECTORS %s float\n',Data.Vectors{k,2});
				fwrite(fid,writeData,'float','b');
			end
		end
	end
	if isfield(Data,'Tensors')
		for k=1:size(Data.Tensors,1)
			if not(isempty(Data.Tensors{k,1}))
				writeData=transpose(Data.Tensors{k,1});
				fprintf(fid,'\nVECTORS %s float\n',Data.Tensors{k,2});
				fwrite(fid,writeData,'float','b');
			end
		end
	end
	fclose(fid);
endfunction