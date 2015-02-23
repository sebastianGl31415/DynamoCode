function [ClebschGordan,IndexOut]=convertWigner2ClebschGordan(Wigner3j,IndexIn)
% convertWigner2ClebschGordan -- Convert Wigner-3j-symbols to 
% Clebsch-Gordan-coefficients. Used to evaluate the K-integral.
%
%     [ClebschGordan,IndexOut]=convertWigner2ClebschGordan(Wigner3j,IndexIn)
%       returns the Clebsch-Gordan-coefficients and related indices for a given 
%       set of Wigner-3j-symbols and corresponding indices.
%
% references -- Arfken, George B. "Mathematical Methods for Physicists"
%            -- http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html

	%% input check
	if not(nargin==2)
		error('2 input arguments are required.')
	elseif not(size(Wigner3j,1)==size(IndexIn,1))
		error('input arguments must have the same length.')
	end
	%% getting index vectors for conversion
	j1=IndexIn(:,1);
	j2=IndexIn(:,2);
	j3=IndexIn(:,3);
	m=IndexIn(:,6);
	%% conversion
	ClebschGordan=(-1).^(j1-j2-m).*sqrt(2*j3+1).*Wigner3j;
	%% index modification
	IndexOut=IndexIn;
	IndexOut(:,6)=-m;
end
