function [Wigner,Ind]=getWignerZeroM(j3,varargin)
% getWignerZeroM -- Calculates all Wigner-3j-Symbols for given j3 . The 
% function returns all values and indices [j1 j2 j3] for which the 
% Wigner-3j-Symbols are not vanishing. The indices for j2 must be stored 
% in the global structure velocity. The computation is performed using 
% recurrence relations for the Wigner-3j-symbols.
%
% The Wigner-3j-Symbols are denoted by standard notation
%
%   / j1  j2  j3 \
%   |            |  and selection rules are applied prior to computing.
%   \ 0   0   0  /
%
%             Wigner = getAllWigner(j3) returns the values as a row vector.
%     [Wigner,Index] = getAllWigner(j3) returns the values as a row vector
%                      and the indices as N-by-6 matrix. Whereas the index 
%                      matrix rows consist of [j1 j2 j3 m1 m2 m3].
%
% If optional 2nd argument is specified as '+' the Wigner-3j-Symbols are 
% computed for j3 + 1. Hence the output is 
%
%   / j1  j2  (j3+1) \
%   |                | .
%   \  0   0     0   /
%
%   See also getWignerZeroM , WignerRecurrence .
	global params
	global velocity
	j2=velocity.l;
	nmax=params.nmax;
	if j3>nmax
		error('1st input (j3) is larger than nmax.');
		return;
	end
	if nargin==1
		[J1,J2,J3,M1,M2,M3]=ndgrid([0:nmax],j2,j3,0,0,0);
	elseif nargin==2
		if strcmp(varargin{1},'+')
			[J1,J2,J3,M1,M2,M3]=ndgrid([0:nmax],j2+1,j3,0,0,0);
		else
			error('Optional 2nd input argument must be the string + .')
		end
	else
		error('Wrong number of input arguments.')
	end
	%% apply selection rule I-1
	Ind=find(bsxfun(@or,bsxfun(@lt,M1,-J1),bsxfun(@gt,M1,J1)));
	%% apply selection rule I-2
	Ind=[Ind; find(bsxfun(@or,bsxfun(@lt,M2,-J2),bsxfun(@gt,M2,J2)))];
	%% apply selection rule I-3
	Ind=[Ind; find(bsxfun(@or,bsxfun(@lt,-M3,-J3),bsxfun(@gt,-M3,J3)))];
	%% apply selection rule II
	Ind=[Ind; find(bsxfun(@ne,M1+M2,-M3))];
	%% apply selection rule III
	Ind=[Ind; find(bsxfun(@or,bsxfun(@gt,J3,J1+J2),bsxfun(@lt,J3,abs(J1-J2))))];
	%% set elements to zero
	[J1,J2,J3,M1,M2,M3]=setInd2Zero(Ind,J1,J2,J3,M1,M2,M3);
	%% find all remaining indices
	orJ=bsxfun(@or,bsxfun(@ne,J1,0), ...
				bsxfun(@or,bsxfun(@ne,J2,0),bsxfun(@ne,J3,0)));
	orM=bsxfun(@or,bsxfun(@ne,M1,0), ...
				bsxfun(@or,bsxfun(@ne,M2,0),bsxfun(@ne,M3,0)));
	Ind=find(bsxfun(@or,orJ,orM));
	%% create vectors of remaining indices
	J1=J1(Ind);
	J2=J2(Ind);
	J3=J3(Ind);
	M1=M1(Ind);
	M2=M2(Ind);
	M3=M3(Ind);
	%% create unique input to recurrence i.e. no doubling of rows
	RecurrenceInput=unique([J2 J3 M1 M2 M3],'rows');
	%% allocation of output
	Wigner=zeros(length(J1),1);
	Ind=zeros(length(J1),6);
	%% row counter l
	l=0;
	%% start recurrence over all rows of RecurrenceInput
	for k=1:size(RecurrenceInput,1)
		%% perfom recurrence calculation
		[Out,j1Out]=WignerRecurrence(RecurrenceInput(k,:));
		%% delete values above nmax
		j1Out=j1Out(j1Out<=nmax);
		Out=Out(j1Out<=nmax);
		%% write results to output i.e. Wigner and Index fields
		Wigner(l+1:l+length(Out))=Out;
		Ind(l+1:l+length(Out),1)=j1Out;
		Ind(l+1:l+length(Out),2:end)=repmat(RecurrenceInput(k,:),length(Out),1);
		%% increase row counter
		l=l+length(Out);
	end
	if nargin==2
		if strcmp(varargin{1},'+')
			Ind(:,2)=Ind(:,2)-1;
		end
	end
end
