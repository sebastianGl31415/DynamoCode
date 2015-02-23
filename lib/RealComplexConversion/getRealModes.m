function [RealModes,IndRealModes]=getRealModes(SinModes,CosModes)
% getSinCosModes -- converts real mode vector two cell arrays SinModes and
% CosModes with following row content
%    Sin/CosModes{row,1}(:) = a^k_l(:) (row vector),
%    Sin/CosModes{row,2} = k           (scalar),
%    Sin/CosModes{row,3} = l           (scalar),
%
% to real mode vector of the form
%     _          _
%    |  a^( 0)_0  |
%    |  b^( 0)_0  |
%    |  a^( 0)_1  |
%    |  b^( 0)_1  |
%    |  a^( 1)_1  |   and an index vector containing (k,l) pair corresponding 
%    |  b^( 1)_1  |   to real mode vector.
%    |      .     |
%    |      .     |
%    |      .     |
%    |  a^( n)_n  |
%    |_ b^( n)_n _|
%

	%% loading global parameter structure
	global params
	nmax=params.nmax;
	nr=params.nr;
	%% number of modes in input vector
	nmodes=(nmax+1)*(nmax+2);
	%% number of modes related to sin
	nsin=0.5*(nmax+1)*nmax;
	%% number of modes related to cos
	ncos=0.5*(nmax+1)*(nmax+2);
	%% checking of input is correct
	if not(size(SinModes,1)==nsin)
		error('Something wrong with the length of 1st input array.')
	end
	if not(size(CosModes,1)==ncos)
		error('Something wrong with the length of 2nd input array.')
	end
	if not(length(SinModes{1,1})==length(CosModes{1,1}))
		error('1st und 2nd input array are not based on the same grid.')
	end
	%% get number of grid points
	np=length(SinModes{1,1});
	%% array allocation
	RealModes=zeros(nmodes*nr,1);
	IndRealModes=zeros(nmodes,2);
	%% fill RealModes / every 2nd mode is a cos modes
	for k=1:ncos
		%% formula for corresponding row, row = l*(l+1)+k*shift+1
		%% shift is equal 2 to get every second row
		row=CosModes{k,3}*(CosModes{k,3}+1)+CosModes{k,2}*2+1;
% 		if row==3 || row==7 || row==13
% 			keyboard
% 		end
		%% store values according to scheme
		RealModes((row-1)*np+1:row*np)=CosModes{k,1};
		IndRealModes(row,1)=CosModes{k,2};
		IndRealModes(row,2)=CosModes{k,3};
	end
	%% fill RealModes / every 2nd mode is a sin modes accept SinModes{k,3}==0
	%% kick out SinModes{k,2}==0 cases to avoid if in loop
	kickOut=find(cell2mat(SinModes(:,2))==0);
	if not(isempty(kickOut))
		keepIn=find(not(cell2mat(SinModes(:,2))==0));
		SinModes=SinModes(keepIn,:);
		nsin=size(SinModes,1);
	end
	%% for loop over remaining sin modes
	for k=1:nsin
		%% formula for corresponding row, row = l*(l+1)+(k+1)*shift
		%% shift is equal 2 to get every second row
		row=SinModes{k,3}*(SinModes{k,3}+1)+(SinModes{k,2}+1)*2;
		%% store values according to scheme
		RealModes((row-1)*np+1:row*np)=SinModes{k,1};
		IndRealModes(row,1)=SinModes{k,2};
		IndRealModes(row,2)=SinModes{k,3};
	end
	%% find index values equal zero in second index
	row=find(IndRealModes(3:end,2)==0)+2;
	%% correct the index values in second index
	IndRealModes(row,:)=IndRealModes(row-1,:);
end