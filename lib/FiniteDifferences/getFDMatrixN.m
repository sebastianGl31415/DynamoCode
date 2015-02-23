function [Dout]=getFDMatrixN(n,Deriv,Type)
% getFDMatrixN --

	%% input check
	if not(ischar(Deriv)) || not(ischar(Type))
		error('2nd and 3rd inputs must be a characters.')
	elseif not(strcmp(Deriv,'1') || strcmp(Deriv,'2'))
		error('2nd input must be a characters either 1 or 2.')
	elseif not(strcmp(Type,'ext') || strcmp(Type,'red'))
		error('3rd input must be a characters either m or n.')
	end
	%% load global structure
	global params;
	nr=params.nr;
	h=params.h;
	%% function handle for FD matrix
	switch Deriv
		case '1'
			D=@(n) D1Matrix2ndOrder(n,h);
		case '2'
			D=@(n) D2Matrix2ndOrder(n,h);
		otherwise
			error('something is wrong with 2nd input and switch statement.')
	end
	%% get full matrix (nr+1) x (nr+1)
	Dout=D(nr+1);
	switch Type
		case 'red' % reduced grid
			%% delete last row
			Dout=Dout(1:end-1,:);
		case 'ext' % extended grid
			opts.derivative=Deriv;
			opts.locationRobinBoundary='right';
			opts.funcRobinBoundary=@(r) [1+0*r; -(n+1)+0*r];
			opts.pointRobinBoundary=1.0;
			Dout=applyRobinBoundary(Dout,opts);
		otherwise
		error('something is wrong with 3rd input and switch statement.')
	end
end