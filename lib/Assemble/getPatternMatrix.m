function [Pattern,GeoIntegrals]=getPatternMatrix()
% getPatternMatrix -- Returns coupling pattern and geo integral structure. The 
% global structure params must be initialized to call this script. Especially 
% the velocity modes (velocity.k,velocity.l) and the maximum number of series 
% terms (params.nmax) need to be specified.
%
%    [Pattern,GeoIntegrals]=getPatternMatrix() pattern is a sparse matrix of 
% ones and zeros indicating a coupling of modes. GeoIntegrals is a structure 
% containing KField and LField, which consist of rows containing the following 
% information (linear indexing)
% 
%   (m,n)   (i,j)    (k,l)   integral value
% [  row     col    velmode        K/L      ]
%

	%% loading global structure
	global params
	nmax=params.nmax;
	%% allocation
	Pattern=zeros((nmax+1)^2,(nmax+1)^2);
	KField=zeros((nmax+1)^4,4);
	LField=zeros((nmax+1)^4,4);
	%% row counters
	kRow=0; % for KField
	lRow=0; % for LField
	%% for loops over all (m,n)-pairs
	for n=0:nmax
		for m=-n:n
			%% getting all K-integral values for current pair
			[K,IndK]=getKIntegral(m,n);
			%% if there exists an index combination
			if not(isempty(IndK))
				%% get row in pattern matrix i.e. (m,n) pair
				Row=convertModeSub2ModeInd(m,n);
				%% get columns in pattern matrix i.e. (i,j) pair
				Cols=convertModeSub2ModeInd(IndK(:,4),IndK(:,1));
				%% repeate row value for all column vector entries
				Rows=repmat(Row,size(Cols));
				%% convert velocity mode pair (k,l) to linear index and store in index 
				%% matrix
				VelModes=convertModeSub2ModeInd(IndK(:,5),IndK(:,2));
				%% store all information in KField
				KField(kRow+1:kRow+length(K),:)=[Rows Cols VelModes K];
				%% increase row counter
				kRow=kRow+length(K);
				%% convert row-column-pairs to linear indices
				Ind=sub2ind([(nmax+1)^2,(nmax+1)^2],Rows,Cols);
				%% set value in pattern matrix
				Pattern(Ind)=1.0;
			end
			%% getting all L-integral values for current pair
			[L,IndL]=getLIntegral(m,n);
			if not(isempty(IndL))
				%% get row in pattern matrix i.e. (m,n) pair
				Row=convertModeSub2ModeInd(m,n);
				%% get columns in pattern matrix i.e. (i,j) pair
				Cols=convertModeSub2ModeInd(IndL(:,4),IndL(:,1));
				%% repeate row value for all column vector entries
				Rows=repmat(Row,size(Cols));
				%% convert velocity mode pair (k,l) to linear index and store in index 
				%% matrix
				VelModes=convertModeSub2ModeInd(IndL(:,5),IndL(:,2));
				%% store all information in LField
				LField(lRow+1:lRow+length(L),:)=[Rows Cols VelModes L];
				%% increase row counter
				lRow=lRow+length(L);
				%% convert row-column-pairs to linear indices
				Ind=sub2ind([(nmax+1)^2,(nmax+1)^2],Rows,Cols);
				%% set value in pattern matrix
				Pattern(Ind)=1.0;
			end
		end
	end
	%% apply sparse storage of pattern matrix
	Pattern=sparse(Pattern);
	%% reduce size of fields
	KField=KField(1:kRow,:);
	LField=LField(1:lRow,:);
	%% save fields in GeoIntegrals structure
	GeoIntegrals.KField=KField;
	GeoIntegrals.LField=LField;
end
