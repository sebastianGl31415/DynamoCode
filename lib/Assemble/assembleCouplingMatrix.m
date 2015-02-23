function [CouplingMatrix]=assembleCouplingMatrix()
% assembleCouplingMatrix --
%
	fprintf(['------- assembling coupling matrix -------------------------\n']);
	displayOutput()
	ticid=tic;
	%% get pattern and geo integrals
	ticticid=tic;
	[~,GeoIntegrals]=getPatternMatrix();
	fprintf('  getting geo integrals --- timing: %8.5f\n',toc(ticticid));
	displayOutput()
	%% assembling upper sub matrix
	ticticid=tic;
	[CouplingMatrixM]=assembleM(GeoIntegrals);
	fprintf('  assembling toroidal   --- timing: %8.5f\n',toc(ticticid))
	displayOutput()
	%% assembling lower sub matrix
	ticticid=tic;
	[CouplingMatrixN]=assembleN(GeoIntegrals);
	fprintf('  assembling poloidal   --- timing: %8.5f\n',toc(ticticid))
	displayOutput()
	%% combining both sub matrices
	CouplingMatrix=[CouplingMatrixM; CouplingMatrixN];
	fprintf('  done                  --- timing: %8.5f\n',toc(ticid))
	displayOutput()
end