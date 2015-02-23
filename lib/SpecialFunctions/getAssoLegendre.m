function [P]=getAssoLegendre(Ind,x)

	%% input check
	if Ind(1)>Ind(2)
		P=zeros(size(x));
	elseif Ind(1)<0 || Ind(2)<0
		error('negative index is not allowed.')
	else
		k=Ind(1);
		l=Ind(2);
		%% get all associated legendre polynomials 
		%% multiplied by (-1)^m in octave/MATLAB
		legendrepoly=legendre(l,x);
		legendrepoly=reshape(legendrepoly(k+1,:,:),size(x));
		P=(-1)^k*legendrepoly;
	end
end