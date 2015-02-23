function [Cells]=getCells(dimensions)
% getCells -- returns cell matrix with points number according to scheme
%             below.
% origin has point number 0.
% (r=m*dr , 0              , 0         ) has point number m.
% (r=m*dr , theta=n*dtheta , 0         ) has point number n*nr + m.
% (r=m*dr , theta=n*dtheta , phi=o*dphi) has point number o*(ntheta-1)*(nr-1)
%                                                                + n*(nr-1) + m.
% There is no doubling of points in cartesian coordinates. For the numbering 
% scheme you image a 'spider' web in (x+,z) plane. This web fills only the 
% positve x region of the plane. This web is then rotated above z-axis i.e. 
% phi-direction. 
%
%  [Cells]=getCells(dimensions)
%
% author: Sebastian Glane

	%% input check
	if not(isvector(dimensions)) && not(length(dimensions)==3)
		error('2nd input has wrong dimension. ')
	else
		nr=dimensions(1);
		ntheta=dimensions(2);
		nphi=dimensions(3);
	end
	%% point and cell numbering according to scheme
	getPointNumber=@(lr,ltheta,lphi)(lphi-1)*(ntheta-1)*(nr-1)+...
							(ltheta-1)*(nr-1)+lr-min((lphi-1)*(lphi-2),1)*(lphi-2)*(nr-1);
	%% allocation
	Cells=zeros((ntheta-1)*(nr-1)*nphi,8);
	%% row counter
	row_cnt=0;
	%% cells with r=0 and points theta=0 / theta=pi, tetrahedras
	lr=1;
	for lphi=1:(nphi-1)
		row_cnt=row_cnt+1;
		ltheta=1;
		Cells(row_cnt,1)=getPointNumber(lr,ltheta,1);
		Cells(row_cnt,2)=getPointNumber(lr,ltheta+1,lphi);
		Cells(row_cnt,3)=getPointNumber(lr,ltheta+1,lphi+1);
		Cells(row_cnt,4)=0;	% column ordering due to VTK_TETRA
		ltheta=ntheta-1;
		row_cnt=row_cnt+1;
		Cells(row_cnt,1)=getPointNumber(lr,ltheta,lphi);
		Cells(row_cnt,2)=getPointNumber(lr,ltheta,lphi+1);
		Cells(row_cnt,3)=getPointNumber(lr,ltheta+1,1);
		Cells(row_cnt,4)=0;	% column ordering due to VTK_TETRA
	end
	% close in azimuthal direction 
	lphi=nphi;
	ltheta=1;
	row_cnt=row_cnt+1;
	Cells(row_cnt,1)=getPointNumber(lr,ltheta,1);
	Cells(row_cnt,2)=getPointNumber(lr,ltheta+1,lphi);
	Cells(row_cnt,3)=getPointNumber(lr,ltheta+1,1);
	Cells(row_cnt,4)=0;	% column ordering due to VTK_TETRA
	ltheta=ntheta-1;
	row_cnt=row_cnt+1;
	Cells(row_cnt,1)=getPointNumber(lr,ltheta,lphi);
	Cells(row_cnt,2)=getPointNumber(lr,ltheta,1);
	Cells(row_cnt,3)=getPointNumber(lr,ltheta+1,1);
	Cells(row_cnt,4)=0;	% column ordering due to VTK_TETRA
	%% cells with r=0 and all other points in interior, pyramides
	lr=1;
	for lphi=1:(nphi-1)
		for ltheta=2:ntheta-2
			row_cnt=row_cnt+1;
			Cells(row_cnt,1)=getPointNumber(lr,ltheta,lphi);
			Cells(row_cnt,2)=getPointNumber(lr,ltheta+1,lphi);
			Cells(row_cnt,3)=getPointNumber(lr,ltheta+1,lphi+1);
			Cells(row_cnt,4)=getPointNumber(lr,ltheta,lphi+1);
			Cells(row_cnt,5)=0;	% column ordering due to VTK_PYRAMID
		end
	end
	% close in azimuthal direction
	lr=1;
	lphi=nphi;
	for ltheta=2:ntheta-2
		row_cnt=row_cnt+1;
		Cells(row_cnt,1)=getPointNumber(lr,ltheta,lphi);
		Cells(row_cnt,2)=getPointNumber(lr,ltheta+1,lphi);
		Cells(row_cnt,3)=getPointNumber(lr,ltheta+1,1);
		Cells(row_cnt,4)=getPointNumber(lr,ltheta,1);
		Cells(row_cnt,5)=0;	% column ordering due to VTK_PYRAMID
	end
	%% cells with points theta=0 or theta=pi, wedges
	for lphi=1:(nphi-1)
		for lr=1:(nr-2)
			ltheta=1;
			row_cnt=row_cnt+1;	% column ordering due to VTK_WEDGE
			Cells(row_cnt,1)=getPointNumber(lr,ltheta+1,lphi);
			Cells(row_cnt,2)=getPointNumber(lr,ltheta+1,lphi+1);
			Cells(row_cnt,3)=getPointNumber(lr,ltheta,1);
			Cells(row_cnt,4)=getPointNumber(lr+1,ltheta+1,lphi);
			Cells(row_cnt,5)=getPointNumber(lr+1,ltheta+1,lphi+1);
			Cells(row_cnt,6)=getPointNumber(lr+1,ltheta,1);
			ltheta=ntheta-1; 
			row_cnt=row_cnt+1;	% column ordering due to VTK_WEDGE
			Cells(row_cnt,1)=getPointNumber(lr,ltheta,lphi);
			Cells(row_cnt,2)=getPointNumber(lr,ltheta,lphi+1);
			Cells(row_cnt,3)=getPointNumber(lr,ltheta+1,1);
			Cells(row_cnt,4)=getPointNumber(lr+1,ltheta,lphi);
			Cells(row_cnt,5)=getPointNumber(lr+1,ltheta,lphi+1);
			Cells(row_cnt,6)=getPointNumber(lr+1,ltheta+1,1);
		end
	end
	% close in azimuth direction
	lphi=nphi;
	for lr=1:(nr-2)
		ltheta=1;
		row_cnt=row_cnt+1;	% column ordering due to VTK_WEDGE
		Cells(row_cnt,1)=getPointNumber(lr,ltheta+1,lphi);
		Cells(row_cnt,2)=getPointNumber(lr,ltheta+1,1);
		Cells(row_cnt,3)=getPointNumber(lr,ltheta,1);
		Cells(row_cnt,4)=getPointNumber(lr+1,ltheta+1,lphi);
		Cells(row_cnt,5)=getPointNumber(lr+1,ltheta+1,1);
		Cells(row_cnt,6)=getPointNumber(lr+1,ltheta,1);
		ltheta=ntheta-1;
		row_cnt=row_cnt+1;	% column ordering due to VTK_WEDGE
		Cells(row_cnt,1)=getPointNumber(lr,ltheta,lphi);
		Cells(row_cnt,2)=getPointNumber(lr,ltheta,1);
		Cells(row_cnt,3)=getPointNumber(lr,ltheta+1,1);
		Cells(row_cnt,4)=getPointNumber(lr+1,ltheta,lphi);
		Cells(row_cnt,5)=getPointNumber(lr+1,ltheta,1);
		Cells(row_cnt,6)=getPointNumber(lr+1,ltheta+1,1);
	end
	%% cells in the interior, hexahedra
	for lphi=1:nphi-1
		for ltheta=2:ntheta-2
			for lr=1:(nr-2)
				row_cnt=row_cnt+1;
				Cells(row_cnt,1)=getPointNumber(lr,ltheta,lphi);
				Cells(row_cnt,2)=getPointNumber(lr+1,ltheta,lphi);
				Cells(row_cnt,3)=getPointNumber(lr+1,ltheta+1,lphi);
				Cells(row_cnt,4)=getPointNumber(lr,ltheta+1,lphi);
				Cells(row_cnt,5)=getPointNumber(lr,ltheta,lphi+1);
				Cells(row_cnt,6)=getPointNumber(lr+1,ltheta,lphi+1);
				Cells(row_cnt,7)=getPointNumber(lr+1,ltheta+1,lphi+1);
				Cells(row_cnt,8)=getPointNumber(lr,ltheta+1,lphi+1);
			end
		end
	end
	% close in azimuthal direction
	lphi=nphi;
	for ltheta=2:ntheta-2
		for lr=1:(nr-2)
			row_cnt=row_cnt+1;
			Cells(row_cnt,1)=getPointNumber(lr,ltheta,lphi);
			Cells(row_cnt,2)=getPointNumber(lr+1,ltheta,lphi);
			Cells(row_cnt,3)=getPointNumber(lr+1,ltheta+1,lphi);
			Cells(row_cnt,4)=getPointNumber(lr,ltheta+1,lphi);
			Cells(row_cnt,5)=getPointNumber(lr,ltheta,1);
			Cells(row_cnt,6)=getPointNumber(lr+1,ltheta,1);
			Cells(row_cnt,7)=getPointNumber(lr+1,ltheta+1,1);
			Cells(row_cnt,8)=getPointNumber(lr,ltheta+1,1);
		end
	end
end
