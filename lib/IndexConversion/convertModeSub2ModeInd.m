function [Ind]=convertModeSub2ModeInd(varargin)
	if nargin==1
		N=varargin{1};
		if (((ne(size(N,1),1) || ne(size(N,2),1)) || ne(N-floor(N),0)) ...
					|| (sign(N)==-1))
			error('If called with one argument, 1st argument must be scalar positive integer.');
			return;
		end
		Ind=zeros((N+1)^2,1);
		k=1;
		Ind(k)=1;
		for n=1:N
			for m=-n:n
				k=k+1;
				Ind(k)=(n)^2+(n+m+1);
			end
		end
	elseif nargin==2
		m=varargin{1};
		n=varargin{2};
		if not(isempty(find(m>n)))
			error('If called with two argument, 1st argument must be less equal 2nd argument.');
			return;
		end
		Ind=n.^2+n+m+1;
	end
end
