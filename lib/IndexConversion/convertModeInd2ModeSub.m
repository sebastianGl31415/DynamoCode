function [m,n]=convertModeInd2ModeSub(Ind)
	if not(isempty(find(Ind<1)))
		error('Input argument must be larger equal 1.');
		return;
	end
	n=floor(sqrt(Ind-1));
	m=Ind-n.^2-n-1;
end
