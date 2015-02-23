function N=getNormalizationConstant(varargin)
	if nargin==2
		m=varargin{1};
		n=varargin{2};
		if length(m)==length(n)
			N=sqrt((2*n+1)/(4*pi).*factorial(n-abs(m))./factorial(n+abs(m)));
		else
			error('Other cases are not implemented yet.')
		end
	else
		error('Other cases are not implemented yet.')
	end
end