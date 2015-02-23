function [RSquaredMatrix]=assembleRSquaredMatrix()
% assembleRSquaredMatrix --
	global params;
	fprintf(['------- assembling r^2 matrix -----------------------------\n']);
	displayOutput()
	tic;
	%% load global structure
	nmax=params.nmax;
	Rhm=params.Rhm;
	Rhn=params.Rhn;
	%% number of modes
	nmodes=(nmax+1)^2;
	%% identity for kronecker products
	Imode=speye(nmodes,nmodes);
	%% creating upper and lower part
	RSquaredM=kron(Imode,Rhm.^2);
	RSquaredN=kron(Imode,Rhn.^2);
	%% joining using block diagonal
	RSquaredMatrix=blkdiag(RSquaredM,RSquaredN);
	fprintf('  done                  --- timing: %8.5f\n',toc)
	displayOutput()
end
