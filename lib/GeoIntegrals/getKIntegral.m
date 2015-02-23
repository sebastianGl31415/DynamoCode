function [K,Index]=getKIntegral(m,n)
% getKIntegral -- Calculates the Adam Gaunt Integral for given m and n. 
% The function returns all values and indices [j l n i k m] for which the 
% Adam Gaunt Integral is not vanishing. The values for l must stored in the 
% global structure params. The Adam Gaunt Integral indexing is
%    ikm
%   K    .
%    jln
%
%	        K = getKIntegral(m,n) returns the values as a row vector.
%	[K,Index] = getKIntegral(m,n) returns the values as a row vector
%             and the indices as N-by-6 matrix. Whereas the index matrix rows 
%             consist of [j l n i k m].
%
%   See also getLIntegral.

	%% get 3j-symbols (standard case)
	[Wigner3j,Index]=getAllWigner(n,m);
	%% get 3j-symbols (zero m case)
	[WignerZeroM,IndexZeroM]=getWignerZeroM(n);
	%% allocation
	Wigner3jZeroM=zeros(size(Wigner3j));
	%% reducing field since last three columns are always zero
	IndexZeroM=IndexZeroM(:,1:3);
	%% find all cases where IndexZeroM and Index are equal
	for k=1:size(IndexZeroM,1)
		[Row,~]=find(all(bsxfun(@eq,Index(:,1:3),IndexZeroM(k,:)),2));
		%% store matching cases in new field
		Wigner3jZeroM(Row)=WignerZeroM(k);
	end
	%% vectorization for evaluation of factor
	j1=Index(:,1);
	j2=Index(:,2);
	j3=Index(:,3);
	m1=Index(:,4);
	m2=Index(:,5);
	m3=Index(:,6);
	%% calculate factor for integral
	Factor=sqrt((2*j1+1).*(2*j2+1)./(4*pi*(2*j3+1)));
	%% convert 3j-symbols to Clebsch-Gordan coefficients (standard case)
	[ClebschGordan,Index]=convertWigner2ClebschGordan(Wigner3j,Index);
	%% convert 3j-symbols to Clebsch-Gordan coefficients (zero m case)
	IndexZeroM=[Index(:,1:3)	zeros(size(Index(:,1:3)))];
	[ClebschGordanZeroM,~]=convertWigner2ClebschGordan(Wigner3jZeroM,IndexZeroM);
	%% calculate K integral values
	K=Factor.*ClebschGordanZeroM.*ClebschGordan;
	%% find non-zero entries
	[Row,~]=find(ne(K,0.0));
	%% delete zero entries
	K=K(Row);
	Index=Index(Row,:);
end
