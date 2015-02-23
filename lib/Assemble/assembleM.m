function [CouplingM]=assembleM(GeoIntegrals)
%	assembleM -- assembles the right-hand sided of the Bullard-Gellman 
% equations, i.e. the coupling matrix for toroidal mode (m).

	%% load global structure
	global params;
	nmax=params.nmax;
	nr=params.nr;
	rhm=params.rhm;
	Rhm=params.Rhm;
	%% finite difference matrices
	D1m=params.D1m;
	Im=speye(nr,nr); % identity matrix for toroidal modes
	In=speye(nr,nr+1); % identity matrix for toroidal modes
	%% submatrices are tridiagonal, hence number of elements in 
	%% nr-by-nr tridiagonal matrix (m case)
	subMatrixElementsM=3*nr-2;
	%% nr-by-(nr+1) tridiagonal matrix (n case)
	subMatrixElementsN=3*nr-1;
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
		MatrixRowStart=convertModeInd2MatrixSub(RowMode,'m');
		%--< p x m , m >-------------------------------------------------------%
		%% get velocity field evaluated on grid
		P=getP(k,l,Rhm);
		%% attention A is multiplied with the m mode
		if not(all(full(diag(P))==0))
			writeA=logical(0);
			if not(n==0) % no division by zero
				A=-K/2*(l*(l+1)/n/(n+1)*getQ(j,n,l)*Rhm*(P*D1m ... 
					+sparse([1:nr],[1:nr],full(D1m*diag(P)),nr,nr)*Im) ...
					+getQ(j,l,n)*Rhm*DOperator(P,'1')*Im);
				writeA=logical(1);
			elseif ((n==0) && (getQ(j,n,l)==0)) % avoid division by zero warning
				A=-K/2*getQ(j,l,n)*Rhm*DOperator(P,'1')*Im;
				writeA=logical(1);
			else % do not know if this ever happens
				error('There is a real division by zero.')
			end
			if writeA
				%% get corresponding matrix columns for m
				MatrixColStart=convertModeInd2MatrixSub(ColMode,'m');
				%% rows columns of submatrix and values
				[subRow,subCol,subVal]=find(A);
				%% setting row column values vector of large matrix
				nrows=length(subRow);
				rowK(row_cnt+1:row_cnt+nrows)=MatrixRowStart+subRow-1;
				colK(row_cnt+1:row_cnt+nrows)=MatrixColStart+subCol-1;
				valK(row_cnt+1:row_cnt+nrows)=subVal;
				%% increase row counter
				row_cnt=row_cnt+nrows;
			end
		end
		%--< o x n , m >-------------------------------------------------------%
		%% get velocity field evaluated on grid
		O=getO(k,l,Rhm);
		if not(all(full(diag(O))==0))
			writeB=logical(0);
			%% get correct finite difference matrix for poloidal modes
			D1n=params.D1n(n,'red');
			%% attention B is multiplied with poloidal mode n
			if not(n==0) % no division by zero
				B=K/2*(j*(j+1)/n/(n+1)*getQ(l,n,j)*Rhm*(O*D1n ...
					+sparse([1:nr],[1:nr],full(D1m*diag(O)),nr,nr)*In) ... 
					+getQ(j,l,n)*O*(Rhm*D1n+In));
				writeB=logical(1);
			elseif ((n==0) && (getQ(l,n,j)==0)) % avoid division by zero warning
				B=K/2*getQ(j,l,n)*O*(Rhm*D1n+In);
				writeB=logical(1);
			else % do not know if this ever happens
				error('There is a real division by zero.')
			end
			if writeB
				%% get corresponding matrix columns for n
				MatrixColStart=convertModeInd2MatrixSub(ColMode,'n');
				%% rows columns of submatrix and values
				[subRow,subCol,subVal]=find(B);
				%% setting row column values vector of large matrix
				nrows=length(subRow);
				rowK(row_cnt+1:row_cnt+nrows)=MatrixRowStart+subRow-1;
				colK(row_cnt+1:row_cnt+nrows)=MatrixColStart+subCol-1;
				valK(row_cnt+1:row_cnt+nrows)=subVal;
				%% increase row counter
				row_cnt=row_cnt+nrows;
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
		MatrixRowStart=convertModeInd2MatrixSub(RowMode,'m');
		%%--< o x m , m >-------------------------------------------------------%
		%% get velocity field evaluated on grid
		O=getO(k,l,Rhm);
		if not(all(full(diag(O))==0))
			%% attention A is multiplied with the m mode
			A=-L*Rhm*O*Im;
			%% get corresponding matrix columns for m
			MatrixColStart=convertModeInd2MatrixSub(ColMode,'m');
			%% rows columns of submatrix and values
			[subRow,subCol,subVal]=find(A);
			%% setting row column values vector of large matrix
			nrows=length(subRow);
			rowL(row_cnt+1:row_cnt+nrows)=MatrixRowStart+subRow-1;
			colL(row_cnt+1:row_cnt+nrows)=MatrixColStart+subCol-1;
			valL(row_cnt+1:row_cnt+nrows)=subVal;
			%% increase row counter
			row_cnt=row_cnt+nrows;
		end
		%%--< p x n , m >-------------------------------------------------------%
		%% get velocity field evaluated on grid
		P=getP(k,l,Rhm);
		if not(all(full(diag(P))==0))
			writeB=logical(0);
			%% get correct finite difference matrices for poloidal modes
			D1n=params.D1n(n,'red');
			D2n=params.D2n(n,'red');
			%% attention B is multiplied with the n mode
			if not(n==0)
				B=L*(-DOperator(P,'1')*(Rhm*D1n+In) ... 
					+l*(l+1)/n/(n+1)*( ...
						DOperator(diag(P)./rhm,'1')*(Rhm.^2*D1n+Rhm*In) ... 
							+P*(Rhm*D2n+2*D1n)) ...
					+j*(j+1)/n/(n+1)*Rhm*(DOperator(P,'1')*D1n ...
						+DOperator(P,'2')*In));
				writeB=logical(1);
			else
				error('There is a division by zero.')
			end
			if writeB
				%% get corresponding matrix columns for n
				MatrixColStart=convertModeInd2MatrixSub(ColMode,'n');
				%% rows columns of submatrix and values
				[subRow,subCol,subVal]=find(B);
				%% setting row column values vector of large matrix
				nrows=length(subRow);
				rowL(row_cnt+1:row_cnt+nrows)=MatrixRowStart+subRow-1;
				colL(row_cnt+1:row_cnt+nrows)=MatrixColStart+subCol-1;
				valL(row_cnt+1:row_cnt+nrows)=subVal;
				%% increase row counter
				row_cnt=row_cnt+nrows;
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
	CouplingM=sparse([rowK rowL],[colK colL],[valK valL],...
											(nmax+1)^2*nr,(nmax+1)^2*(2*nr+1));
end