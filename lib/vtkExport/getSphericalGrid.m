function [R,Theta,Phi]=getSphericalGrid(radius,dimensions)
% getSphericalGrid -- returns (R,Theta,Phi) vectors of spherical coordinates.
%
%  [R,Theta,Phi]=getSphericalGrid(radius,dimensions)
%
% author: Sebastian Glane

%% input check
	if not(isvector(dimensions)) && not(length(dimensions)==3)
		error('2nd input has wrong dimension. ')
	else
		nr=dimensions(1);
		ntheta=dimensions(2);
		nphi=dimensions(3);
	end
	r=linspace(0,radius,nr)';
	theta=linspace(0,pi,ntheta)';
	phi=linspace(0,2*pi,nphi+1)';
	[R,Theta,Phi]=ndgrid(r,theta,phi);
	R=reshape(R,[prod(size(R)) 1]);
	Theta=reshape(Theta,[prod(size(Theta)) 1]);
	Phi=reshape(Phi,[prod(size(Phi)) 1]);
end
