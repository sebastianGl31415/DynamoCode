function [SphericalHarmonic]=CosSphericalHarmonic(Ind,theta,phi)
% CosSphericalHarmonics --
	if not(size(Ind,1)==2) && not(size(Ind,2)==2)
		error('1st input must a vector of length 2.')
	elseif Ind(1)>Ind(2)
		error('1st index is larger than 2nd.')
	elseif not(length(theta)==length(phi))
		error('2nd and 3rd must be of the same size.')
	else 
		m=Ind(1);
		n=Ind(2);
	end
	%% get associated legendre polynomial
	SphericalHarmonic=getAssoLegendre([m n],cos(theta));
	%% spherical harmonic
	SphericalHarmonic=cos(m*phi).*SphericalHarmonic;
end