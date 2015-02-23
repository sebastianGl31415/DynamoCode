function [L,Index]=getLIntegral(m,n)
% getLIntegral -- Calculates the Elsasser Integral for given m and n. 
% The function returns all values and indices [j l n i k m] for which the 
% Elsasser Integral is not vanishing. The values for l must stored in the 
% global structure params. The Elsasser Integral indexing is
%    ikm
%   L    .
%    jln
%
%	        L = getLIntegral(m,n) returns the values as a row vector.
%	[L,Index] = getLIntegral(m,n) returns the values as a row vector
%             and the indices as N-by-6 matrix. Whereas the index matrix rows 
%             consist of [j l n i k m].
%
%   See also getKIntegral.

	%% get 3j-symbols (standard case)
	[Wigner3j,Index]=getAllWigner(n,m);
	%% get 3j-symbols (zero m and j2+1 case)
	[WignerZeroM,IndexZeroM]=getWignerZeroM(n,'+');
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
	m3=Index(:,6);
	%% calculate factor for integral
	Factor=((-1).^m3)*i/2.*sqrt((2*j1+1).*(2*j2+1).*(2*j3+1)/4/pi) ... 
			.*sqrt((j1+j2+j3+2).*(j1+j2-j3+1).*(-j1+j2+j3+1).*(j1-j2+j3));
	%% calculate L integral values
	L=Factor.*Wigner3j.*Wigner3jZeroM;
	%% find non-zero entries
	[Row,~]=find(ne(L,0.0));
	%% delete zero entries
	L=L(Row);
	Index=Index(Row,:);
	%% modify m-index for correction
	Index(:,6)=-Index(:,6);
end