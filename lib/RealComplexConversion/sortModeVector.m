function [sortedModeVector,Ind]=sortModeVector(unsortedModeVector,opt)
% sortModeVector -- modifying the mode vectors according to 2nd input argument.
% 
% If 2nd input argument is 'ext' (extend), sortModeVector extends the input 
% vector
%     _          _       _          _
%    |  c^( 0)_0  |     |  c^( 0)_0  | 
%    |  c^(-1)_1  |     |  c^(-0)_0  |
%    |  c^( 0)_1  |     |  c^( 0)_1  |
%    |  c^( 1)_1  |     |  c^(-0)_1  |
%    |      .     |     |  c^( 1)_1  |
%    |      .     | to  |  c^(-1)_1  |  where as the c^(-0)_k are phantom
%    |      .     |     |      .     |  entries equal to zero. This reordering
%    |_ c^( n)_n _|     |      .     |  is due to conversion from complex to 
%                       |      .     |  real.
%                       |      .     |
%                       |_ c^(-n)_n _|
%
% If 2nd input argument is 'red' (reduced), sortModeVector extends the input 
% vector
%     _          _       _          _
%    |  c^( 0)_0  |     |  c^( 0)_0  | 
%    |  c^(-0)_1  |     |  c^(-1)_1  |
%    |  c^( 0)_1  |     |  c^( 0)_1  |
%    |  c^(-0)_1  |     |  c^( 1)_1  |
%    |  c^( 1)_1  |     |      .     |
%    |  c^(-1)_1  | to  |      .     |  where as the phantom entires c^(-0)_k 
%    |      .     |     |      .     |  are eliminated. This reordering is due 
%    |      .     |     |_ c^(-n)_n _|  to conversion from complex to real.
%    |      .     |
%    |      .     |
%    |_ c^(-n)_n _|
%
 
	if not(ischar(opt))
		error('2nd input must be character array.')
	end
	%% load global parameter structure
	global params
	nmax=params.nmax;
	nr=params.nr;
	%% cell array to store mode permutation matrices 
	M=cell(nmax,1);
	if strcmp(opt,'ext')
		nmodes=(nmax+1)^2;
	elseif strcmp(opt,'red')
		nmodes=(nmax+1)*(nmax+2);
	end
	%% allocation of indices k and l
	indk=zeros(nmodes,1);
	indl=zeros(nmodes,1);
	%% check if poloidal or toroidal field by number of grid points
	if (size(unsortedModeVector,1)/nmodes)==nr
		np=nr;
	elseif (size(unsortedModeVector,1)/nmodes)==(nr+1)
		np=nr+1;
	else
		error('Something is wrong with the number of points.')
	end
	%% identity for grid points
	Ir=speye(np);
	if strcmp(opt,'ext')
		%% row counter for k and l
		rowcount=0;
		%% loop over all modes c^k_l for fixed l
		for l=0:nmax
			%% number of submodes
			nsubmodes=2*l+1;
			%% setting k and l values
			indk(rowcount+1:rowcount+nsubmodes)=[-l:l]';
			indl(rowcount+1:rowcount+nsubmodes)=l;
			rowcount=rowcount+nsubmodes;
			%% extended identity matrix, zeros in last row
			Imod=speye(nsubmodes+1,nsubmodes);
			%% permutation vector
			p=getPermutationVector(l+1);
			%% store mode permutation matrix
			M{l+1,1}=Imod(p,:);
		end
	elseif strcmp(opt,'red')
		%% row counter for k and l
		rowcount=0;
		%% loop over all modes c^k_l for fixed l
		for l=0:nmax
			%% number of submodes
			nsubmodes=2*l+2;
			%% setting k and l values
			indk(rowcount+1:rowcount+nsubmodes)=repelems(0:l,[1:l+1; 2*ones(1,l+1)]);
			indk(rowcount+2:2:rowcount+nsubmodes)=-indk(rowcount+2:2:rowcount...
																															+nsubmodes);
			indl(rowcount+1:rowcount+nsubmodes)=l;
			rowcount=rowcount+nsubmodes;
			%% extended identity matrix, zeros in last row
			Imod=speye(nsubmodes-1,nsubmodes);
			%% permutation vector
			p=getPermutationVector(l+1);
			%% store mode permutation matrix
			M{l+1,1}=Imod(:,p);
		end
	end
	%% create overall permutation matrix as block diagonal
	P=blkdiag(M{:});
	%% permutation with kronecker product to extend permutation on grid
	sortedModeVector=kron(P,Ir)*unsortedModeVector;
	%% permutation of indices k and l
	permindk=P*indk;
	permindl=P*indl;
	if strcmp(opt,'ext')
		%% find zeros in permindl occuring through matrix multiply
		row_permindl_zero=find(permindl==0);
		%% correct zero entries
		permindl(row_permindl_zero(2:end))=permindl(row_permindl_zero(2:end)-1);
	end
	%% write output
	Ind=[permindk permindl];
end