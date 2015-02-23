function []=initializeDiscretization()
	global params;
	fprintf('------- discretization information -------------------------\n');
	displayOutput()
	nr=params.nr;
	rh=linspace(0,1,nr+2)';
	params.h=rh(2)-rh(1);
	params.rhm=rh(2:end-1);
	params.rhn=rh(2:end);
	params.Rhm=sparse([1:nr],[1:nr],params.rhm,nr,nr);
	params.Rhn=sparse([1:(nr+1)],[1:(nr+1)],params.rhn,(nr+1),(nr+1));
	params.D1hm=D1Matrix2ndOrder(nr,params.h);
	params.D2hm=D2Matrix2ndOrder(nr,params.h);
	params.D1hn=D1Matrix2ndOrder(nr+1,params.h);
	params.D2hn=D2Matrix2ndOrder(nr+1,params.h);
	params.D1m=params.D1hm;
	params.D2m=params.D2hm;
	params.D1n=@(n,gridType) getFDMatrixN(n,'1',gridType);
	params.D2n=@(n,gridType) getFDMatrixN(n,'2',gridType);
	fprintf(['  number of points (nr)  -- %d \n',...
'  number of modes (nmax) -- %d\n',...
'  total number of modes  -- %d\n',...
'  grid spacing (h)       -- %10.9f\n',...
'  finite diffence scheme -- %s\n',...
'  problem size           -- %d\n'],nr,params.nmax,...
		2*(params.nmax+1)^2,params.h,params.diffMethod,...
		(params.nmax+1)^2*(2*nr+1));
	displayOutput()
end
