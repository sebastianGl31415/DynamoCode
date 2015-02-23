function [Ind]=getModeIndexComplex(m,n)
	global params;
	nr=params.nr;
	Row=convertModeSub2ModeInd(m,n);
	Ind=[Row*nr+1:(Row+1)*nr]';
end