function [SinModes,CosModes]=getSinCosModes(RealModes,IndRealModes)
% getSinCosModes -- converts real mode vector
%     _          _  
%    |  a^( 0)_0  | to two cell arrays SinModes and CosModes with following
%    |  b^( 0)_0  | row content
%    |  a^( 0)_1  |
%    |  b^( 0)_1  |      Sin/CosModes{row,1}(:) = a^k_l(:) (row vector),
%    |  a^( 1)_1  |      Sin/CosModes{row,2} = k           (scalar),
%    |  b^( 1)_1  |      Sin/CosModes{row,3} = l           (scalar).
%    |      .     |
%    |      .     | Index vector for input needs to be specified as 2nd input.
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
	%% check if poloidal or toroidal field by number of grid points
	if (size(RealModes,1)/nmodes)==nr
		%% identity for grid points
		np=nr;
	elseif (size(RealModes,1)/nmodes)==(nr+1)
		%% identity for grid points
		np=(nr+1);
	else
		error('Something is wrong with the number of points.')
	end
	%% array allocation
	SinModes=cell(nsin,3);
	CosModes=cell(ncos,3);
	%% fill CosModes / every 2nd mode is a cos modes
	for k=1:2:nmodes
		CosModes{(k+1)/2,1}=RealModes((k-1)*np+1:k*np);
		CosModes{(k+1)/2,2}=IndRealModes(k,1);
		CosModes{(k+1)/2,3}=IndRealModes(k,2);
	end
	%% find the rows where sin is zero / a^0_l and b^0_l are eliminated
	sin_rows=find(not(IndRealModes(:,1)==0));
	%% every 2nd row is a sin row / sin_rows is [a^k_l b^k_l ... ]
	sin_rows=sin_rows(2:2:end);
	%% loop over remaining rows
	for k=1:length(sin_rows)
		SinModes{k,1}=RealModes((sin_rows(k)-1)*np+1:sin_rows(k)*np);
		SinModes{k,2}=IndRealModes(sin_rows(k),1);
		SinModes{k,3}=IndRealModes(sin_rows(k),2);
	end
end