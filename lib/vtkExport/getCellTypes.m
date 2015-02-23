function [CellsOut,CellTypes]=getCellTypes(Cells)
% getCellTypes -- returns Cells as a cell array and matrix CellTypes. The cell 
%    array Cells consists of the different cell types i.e. tetrahedrons, 
%    pyramides, wedges and hexahedrons. The matirx CellTypes consists fo the 
%    VTK-identifiers for different cell types. This reordering is due to output 
%    format.
%
%  [CellsOut,CellTypes]=getCellTypes(Cells)
%
% author: Sebastian Glane

	%% allocation
	CellsOut=cell(4,1);
	CellTypes=zeros(size(Cells,1),1);
	%% row counter for CellTypes
	row_cnt=0;
	%% find all tetrahedrons
	[row,~]=find(sum(Cells(:,4:8),2)==0);
	CellsOut{1}=[4*ones(size(row)) Cells(row,1:4)];
	Cells(row,:)=[];	% remove tetrahedrons rows
	CellTypes(1:length(row))=10; % vtk id for VTK_TETRA
	row_cnt=row_cnt+length(row);
	%% find all pyramides
	[row,~]=find(sum(Cells(:,5:8),2)==0);
	CellsOut{2}=[5*ones(size(row)) Cells(row,1:5)];
	Cells(row,:)=[];	% remove pyramides rows
	CellTypes(row_cnt+1:row_cnt+length(row))=14; % vtk id for VTK_PYRAMID
	row_cnt=row_cnt+length(row);
	%% find all wedges
	[row,~]=find(sum(Cells(:,7:8),2)==0); % 
	CellsOut{3}=[6*ones(size(row)) Cells(row,1:6)];
	Cells(row,:)=[];	% remove wedges rows
	CellTypes(row_cnt+1:row_cnt+length(row))=13; % vtk id for VTK_WEDGE
	row_cnt=row_cnt+length(row);
	%% rest rows are pyramides
	CellsOut{4}=[8*ones(size(Cells,1),1) Cells];
	CellTypes(row_cnt+1:row_cnt+size(Cells,1))=12; % vtk id for VTK_HEXAHEDRON
end
