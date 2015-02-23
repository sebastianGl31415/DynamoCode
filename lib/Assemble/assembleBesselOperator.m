function [BesselOperator]=assembleBesselOperator()
% assembleBesselOperator --
	global params;
	fprintf([	'------- assembling bessel operator matrix------------------\n']);
	displayOutput()
	tic;
	%% load global structure
	nmax=params.nmax;
	nr=params.nr;
	Rhm=params.Rhm;
	Rhn=params.Rhn;
	%% identity matrices for m and n case
	Im=speye(nr,nr);
	In=speye(nr+1,nr+1);
	%% load finite difference matrix for m
	D1m=params.D1m;
	D2m=params.D2m;
	%% number of modes
	nmodes=(nmax+1)^2;
	%% sparse matrix with n on diagonal corresponding to mode
	n=floor(sqrt([1:nmodes]'-1));
	N=sparse([1:nmodes],[1:nmodes],n,nmodes,nmodes);
	%% identity for kronecker products
	Imode=speye(nmodes,nmodes);
	%%-- assembling upper sub matrix case (m) ----------------------------------%
	BesselOperatorM=kron(Imode,Rhm.^2*D2m+2*Rhm*D1m)-kron(N.*(N+1),Im);
	if not(issparse(BesselOperatorM))
		BesselOperatorM=sparse(BesselOperatorM);
	end
	%%-- assembling lower sub matrix case (n) ----------------------------------%
	%% submatrices are tridiagonal, hence number of elements in 
	%% (nr+1)-by-(nr+1) tridiagonal matrix (n case)
	subMatrixElements=3*(nr+1)-2;
	%% maximum number of elements for KField corresponding submatrices
	numelBesselN=nmodes*subMatrixElements;
	%% allocation of vectors for rows columns and values
	rowBesselN=zeros(1,numelBesselN);
	colBesselN=zeros(1,numelBesselN);
	valBesselN=zeros(1,numelBesselN);
	%% row counter
	row_cnt=0;
	for k=1:nmodes
		%% getting FD matrix modified for current 
		D1n=params.D1n(n(k),'ext');
		D2n=params.D2n(n(k),'ext');
		%% Bessel operator submatrix
		A=Rhn.^2*D2n+2*Rhn*D2n-n(k)*(n(k)+1)*In;
		%% rows in large matrix
		MatrixRowStart=convertModeInd2MatrixSub(k,'n');
		%% row correction
		MatrixRowStart=MatrixRowStart-(nmax+1)^2*nr;
		%% rows columns of submatrix and values
		[subRow,subCol,subVal]=find(A);
		%% setting row column values vector of large matrix
		nrows=length(subRow);
		rowBesselN(row_cnt+1:row_cnt+nrows)=MatrixRowStart+subRow-1;
		colBesselN(row_cnt+1:row_cnt+nrows)=MatrixRowStart+subCol-1;
		valBesselN(row_cnt+1:row_cnt+nrows)=subVal;
		%% increase row counter
		row_cnt=row_cnt+nrows;
	end
	%% check if remaining rows are all zeros
	if not(all(rowBesselN(row_cnt+1:end)==0))
		error('error in rowBesselN.')
	elseif not(all(colBesselN(row_cnt+1:end)==0))
		error('error in colBesselN.')
	elseif not(all(valBesselN(row_cnt+1:end)==0))
		error('error in valBesselN.')
	end
	%% cut vectors at lower end
	rowBesselN=rowBesselN(1:row_cnt);
	colBesselN=colBesselN(1:row_cnt);
	valBesselN=valBesselN(1:row_cnt);
	%% large sparse BesselOperatorN matrix
	BesselOperatorN=sparse(rowBesselN,colBesselN,valBesselN,...
										(nmax+1)^2*(nr+1),(nmax+1)^2*(nr+1));
	%%-- assembling BesselOperator as block diagonal matrix --------------------%
	BesselOperator=blkdiag(BesselOperatorM,BesselOperatorN);
	fprintf('  done                  --- timing: %8.5f\n',toc);
	displayOutput();
end
