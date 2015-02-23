function [MatrixIndexStart]=convertModeInd2MatrixSub(ModeIndex,Type)
% convertModeInd2MatrixSub -- converts 
%
	%% input check
	if not(ischar(Type))
		error('2nd input must be a string.')
	elseif not(strcmp(Type,'m') | strcmp(Type,'n'))
		error('2nd input must be a string either m or n.')
	end
	%% load params
	global params
	nr=params.nr;
	nmax=params.nmax;
	%% 
	if strcmp(Type,'m')
		MatrixIndexStart=(ModeIndex-1)*nr+1;
	elseif strcmp(Type,'n')
		MatrixIndexStart=(nmax+1)^2*nr+(ModeIndex-1)*(nr+1)+1;
	end
end