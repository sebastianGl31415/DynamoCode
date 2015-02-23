function [RealModes,IndRealModes]=convertComplex2Real(ComplexModes)
% convertComplex2Real --
	
	global params
	nmax=params.nmax;
	nr=params.nr;
	if not(size(ComplexModes,1)==(nmax+1)^2*nr || ...
			size(ComplexModes,1)==(nmax+1)^2*(nr+1))
		error('The length of the input is not correct.')
	elseif not(size(ComplexModes,2)==1)
		error('Only vector inputs are accepted.')
	end
	%% get sorted mode vector
	[sortedComplexModes,Ind]=sortModeVector(ComplexModes,'ext');
	k=Ind(:,1);
	l=Ind(:,2);
	%% normalization constant
	N=getNormalizationConstant(k,l);
	%% number modes of sorted vector (observe that nmodes is always even)
	nmodes=(nmax+1)*(nmax+2);
	%% sparse diagonal matrix with every 2nd entry N^k_l on diagonal
	%% every 2nd entry since N^k_l = N^(-k)_l
	spdiagN=sparse([1:nmodes/2],[1:nmodes/2],N(1:2:end),nmodes/2,nmodes/2);
	%% check if poloidal or toroidal field by number of grid points
	if (size(sortedComplexModes,1)/nmodes)==nr
		%% identity for grid points
		np=nr;
	elseif (size(sortedComplexModes,1)/nmodes)==(nr+1)
		%% identity for grid points
		np=nr+1;
	else
		error('Something is wrong with the number of points.')
	end
	Igrid=speye(np);
	%% check imaginary unit not be overwritten
	if not(iscomplex(i))
		error('i is not the imaginary unit anymore')
	end
	%% matrix converting complex pair to real pairs
	Mc2r=[1 1; i -i];
	%% build large conversion matrix by kronecker products
	ConversionMatrix=kron(kron(spdiagN,Mc2r),Igrid);
	%% conversion by sparse matrix vector product
	RealModes=ConversionMatrix*sortedComplexModes;
	%% real modes have only positive indices
	IndRealModes=abs(Ind);
	%% correct the rows corresponding to b^0_l
	%% find rows with k=0, hence all a^0_l and b^0_l
	row=find(IndRealModes(:,1)==0);
	%% every 2nd row belongs to b^0_l
	row=row(2:2:end);
	%% loop over to set all b^0_l equal zero
	for k=1:length(row)
		RealModes((row(k)-1)*np+1:row(k)*np)=0.0;
	end
end