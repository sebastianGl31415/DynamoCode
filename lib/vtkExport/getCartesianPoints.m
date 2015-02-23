function [Points,delRows]=getCartesianPoints(R,Theta,Phi)
% getCartesianPoints -- returns cartesian points for given spherical 
%    coordinates. No doubling of points is allowed.
%
%  [Points,delRows]=getCartesianPoints(R,Theta,Phi)
%
% author: Sebastian Glane

	%% ensure that input are row vectors
	if not(size(R,1)==length(R))
		R=R';
	end
	if not(size(Theta,1)==length(Theta))
		Theta=Theta';
	end
	if not(size(Phi,1)==length(Phi))
		Phi=Phi';
	end
	%% spherical coordinates
	X=R.*cos(Phi).*sin(Theta);
	Y=R.*sin(Phi).*sin(Theta);
	Z=R.*cos(Theta);
	%% set of all points
	Points=[X Y Z];
	%% find point components numerically equal to zero
	tol=100*eps;
	Ind=find(abs(Points)<tol);
	%% set point components numerically equal to zero to exactly zero
	Points(Ind)=0.0;
	%% create a unique point set
	[~,Ind]=unique(Points,'rows','first');
	[row,~]=ind2sub(size(Points),sort(Ind));
	delRows=setdiff([1:length(Points)]',row);
	Points=Points(row,:);
end
