function [CartesianField]=SphericalField2Cartesian(Field,SphericalPoints)
	if not(size(SphericalPoints,2)==3)
		error('2nd input must be a n-by-3 matrix containing R,Theta,Phi.')
	end
	%% number of points
	npoints=length(SphericalPoints);
	%% getting z-coordinates
	R=SphericalPoints(:,1);
	Theta=SphericalPoints(:,2);
	Phi=SphericalPoints(:,3);
	%% allocation
	Er=zeros(npoints,3);
	Etheta=zeros(npoints,3);
	Ephi=zeros(npoints,3);
	%% spherical unit vector in cartesian system
	Er(:,1)=cos(Phi).*sin(Theta);
	Er(:,2)=sin(Phi).*sin(Theta);
	Er(:,3)=cos(Theta);
	Etheta(:,1)=cos(Phi).*cos(Theta);
	Etheta(:,2)=sin(Phi).*cos(Theta);
	Etheta(:,3)=-sin(Theta);
	Ephi(:,1)=-sin(Phi);
	Ephi(:,2)=cos(Phi);
	%% extending field physical components
	Field_r=repmat(Field(:,1),[1 3]);
	Field_theta=repmat(Field(:,2),[1 3]);
	Field_phi=repmat(Field(:,3),[1 3]);
	%% cartesian by point wise multiplication
	CartesianField=Field_r.*Er+Field_theta.*Etheta+Field_phi.*Ephi;
end