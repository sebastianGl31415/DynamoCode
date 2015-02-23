function [CouplingN]=assembleN(GeoIntegrals)
% assembleM -- assembles the right-hand sided of the Bullard-Gellman 
% equations, i.e. the coupling matrix for toroidal mode (m).

	%% load global structure
	global params;
	nmax=params.nmax;
	nr=params.nr;
	rhn=params.rhn;
	Rhn=params.Rhn;
	%% finite difference matrices
	D1hn=params.D1hn;
	Im=speye(nr+1,nr); % identity matrix for toroidal modes
	In=speye(nr+1,nr+1); % identity matrix for toroidal modes
	%% submatrices are tridiagonal, hence number of elements in 
	%% nr-by-nr tridiagonal matrix (m case)
	subMatrixElementsM=3*nr-1;
	%% nr-by-(nr+1) tridiagonal matrix (n case)
	subMatrixElementsN=3*(nr+1)-2;
	%% to over estimated we take the maximum
	maxsubMatrixElements=max([subMatrixElementsM subMatrixElementsN]);
	%% getting KField
	KField=GeoIntegrals.KField;
	%% maximum number of elements for KField corresponding submatrices
	numelK=size(KField,1)*maxsubMatrixElements;
	%% allocation of vectors for rows columns and values
	rowK=zeros(1,numelK);
	colK=zeros(1,numelK);
	valK=zeros(1,numelK);
	%% row counter
	row_cnt=0;
	%% for loop over all K-integrals
	for ind=1:size(KField,1)
		%% get mode indices
		RowMode=KField(ind,1);
		ColMode=KField(ind,2);
		%% get (i,j), (k,l) and (m,n) values for formulas below
		[i,j]=convertModeInd2ModeSub(ColMode);
		[k,l]=convertModeInd2ModeSub(KField(ind,3));
		[m,n]=convertModeInd2ModeSub(RowMode);
		%% get K-integral value
		K=KField(ind,4);
		%% rows in coupling matrix
		MatrixRowStart=convertModeInd2MatrixSub(RowMode,'n');
		%% row correction
		MatrixRowStart=MatrixRowStart-(nmax+1)^2*nr;
		%--< p x n , n >-------------------------------------------------------%
		%% get velocity field evaluated on grid
		P=getP(k,l,Rhn);
		if not(all(full(diag(P))==0))
			%% get correct finite difference matrices for poloidal modes
			D1n=params.D1n(n,'ext');
			%% attention A is multiplied with the n mode
			if not(n==0) % no division by zero
				A=K/2/n/(n+1)*(-l*(l+1)*getQ(j,n,l)*P*(Rhn*D1n+In) ...
					+j*(j+1)*getQ(l,n,j)*Rhn*DOperator(P,'1')*In);
				%% get corresponding matrix columns for n
				MatrixColStart=convertModeInd2MatrixSub(ColMode,'n');
				%% rows columns of submatrix and values
				[subRow,subCol,subVal]=find(A);
				%% setting row column values vector of large matrix
				nrows=length(subRow);
				rowK(row_cnt+1:row_cnt+nrows)=MatrixRowStart+subRow-1;
				colK(row_cnt+1:row_cnt+nrows)=MatrixColStart+subCol-1;
				valK(row_cnt+1:row_cnt+nrows)=subVal;
				%% increase row counter
				row_cnt=row_cnt+nrows;
			elseif (((n==0) && (getQ(j,n,l)==0)) && (getQ(l,n,j)==0))	% avoid divi-
				A=0.0;																									% sion by zero
			else % do not know if this ever happens
				error('There is a real division by zero.')
			end
		end
	end
	%% check if remaining rows are all zeros
	if not(all(rowK(row_cnt+1:end)==0))
		error('error in rowK.')
	elseif not(all(colK(row_cnt+1:end)==0))
		error('error in colK.')
	elseif not(all(valK(row_cnt+1:end)==0))
		error('error in valK.')
	end
	%% cut vectors at lower end
	rowK=rowK(1:row_cnt);
	colK=colK(1:row_cnt);
	valK=valK(1:row_cnt);
	%% getting LField
	LField=GeoIntegrals.LField;
	%% maximum number of elements for LField corresponding submatrices
	numelL=size(LField,1)*maxsubMatrixElements;
	%% allocation of vectors for rows columns and values
	rowL=zeros(1,numelK);
	colL=zeros(1,numelK);
	valL=zeros(1,numelK);
	%% reset row counter
	row_cnt=0;
	for ind=1:size(LField,1)
		%% get mode indices
		RowMode=LField(ind,1);
		ColMode=LField(ind,2);
		%% get (i,j), (k,l) and (m,n) values for formulas below
		[i,j]=convertModeInd2ModeSub(ColMode);
		[k,l]=convertModeInd2ModeSub(LField(ind,3));
		[m,n]=convertModeInd2ModeSub(RowMode);
		%% get K-integral value
		L=LField(ind,4);
		%% rows in coupling matrix
		MatrixRowStart=convertModeInd2MatrixSub(RowMode,'n');
		%% row correction
		MatrixRowStart=MatrixRowStart-(nmax+1)^2*nr;
		%%--< o x n , n >-------------------------------------------------------%
		%% get velocity field evaluated on grid
		O=getO(k,l,Rhn);
		if not(all(full(diag(O))==0))
			%% attention A is multiplied with the n mode
			if not(n==0)
				A=-L*j*(j+1)/n/(n+1)*Rhn*O*In;
				%% get corresponding matrix columns for n
				MatrixColStart=convertModeInd2MatrixSub(ColMode,'n');
				%% rows columns of submatrix and values
				[subRow,subCol,subVal]=find(A);
				%% setting row column values vector of large matrix
				nrows=length(subRow);
				rowL(row_cnt+1:row_cnt+nrows)=MatrixRowStart+subRow-1;
				colL(row_cnt+1:row_cnt+nrows)=MatrixColStart+subCol-1;
				valL(row_cnt+1:row_cnt+nrows)=subVal;
				%% increase row counter
				row_cnt=row_cnt+nrows;
			else
				error('There is a division by zero.')
			end
		end
		%%--< p x m , n >-------------------------------------------------------%
		%% get velocity field evaluated on grid
		P=getP(k,l,Rhn);
		if not(all(full(diag(P))==0))
			%% attention B is multiplied with the m mode
			if not(n==0)
				B=-L*l*(l+1)/n/(n+1)*Rhn*P*Im;
				%% get corresponding matrix columns for m
				MatrixColStart=convertModeInd2MatrixSub(ColMode,'m');
				%% rows columns of submatrix and values
				[subRow,subCol,subVal]=find(B);
				%% setting row column values vector of large matrix
				nrows=length(subRow);
				rowL(row_cnt+1:row_cnt+nrows)=MatrixRowStart+subRow-1;
				colL(row_cnt+1:row_cnt+nrows)=MatrixColStart+subCol-1;
				valL(row_cnt+1:row_cnt+nrows)=subVal;
				%% increase row counter
				row_cnt=row_cnt+nrows;
			else
				error('There is a division by zero.')
			end
		end
	end
	%% check if remaining rows are all zeros
	if not(all(rowL(row_cnt+1:end)==0))
		error('error in rowL.')
	elseif not(all(colL(row_cnt+1:end)==0))
		error('error in colL.')
	elseif not(all(valL(row_cnt+1:end)==0))
		error('error in valL.')
	end
	%% cut vectors at lower end
	rowL=rowL(1:row_cnt);
	colL=colL(1:row_cnt);
	valL=valL(1:row_cnt);
	%% create sparse coupling matrix
	CouplingN=sparse([rowK rowL],[colK colL],[valK valL],...
											(nmax+1)^2*(nr+1),(nmax+1)^2*(2*nr+1));
end