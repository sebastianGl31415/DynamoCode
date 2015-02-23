function [du]=DOperator(u,Deriv)
	%% loading global structure
	global params
	nr=params.nr;
	%% check if input is vector or diagonal matrix
	if not(isvector(u)) && not(isscalar(u)) % matrix input
		%% convert diagonal matrix to vector
		u=full(diag(u));
	end
	%% checking grid and getting FD matrices
	if length(u)==nr % toroidal grid
		n=nr;
		rh=params.rhm;
		D1=params.D1hm;
		D2=params.D2hm;
	elseif length(u)==(nr+1) % poloidal grid
		n=nr+1;
		rh=params.rhn;
		D1=params.D1hn;
		D2=params.D2hn;
	end
	%% computing result
	if all(not(u==0))
		switch Deriv
			case '1' % 1st discrete derivative of operator
				du=D1*u+u./rh;
			case '2' % 2nd discrete derivative of operator
				du=D2*u+(D1*u)./rh;
			otherwise
				error('2nd input must be a string 1 or 2.')
		end
	elseif ischar(Deriv) % return zeros if input is zero vector
		du=zeros(size(D1,1),1);
	else
		error('2nd input must be a string 1 or 2.')
	end
	%% create sparse diagonal matrix for output
	du=sparse([1:n],[1:n],du,n,n);
end