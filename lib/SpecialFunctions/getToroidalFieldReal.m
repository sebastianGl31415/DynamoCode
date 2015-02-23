function [Field]=getToroidalFieldReal(Ind,SphericalPoints,varargin)
	%% input check
	if not(size(SphericalPoints,2)==3)
		error('2nd input must be a n-by-3 matrix containing R,Theta,Phi.')
	end
	if nargin==3 && iscell(varargin{1})
		if length(varargin{1})==2
			funcIn=varargin{1};
			for k=1:length(funcIn)
				if not(isa(funcIn{k},'function_handle'))
					error('3rd input must be cell array containing function handles.')
				end
			end
			f_r_cos=funcIn{2};
			f_r_sin=funcIn{2};
		else
			error('3rd input must be a 2-by-1 or 1-by-2 cell array.')
		end
	elseif nargin==3
		error('3rd input must be a cell array.')
	end
	%% number of points
	npoints=length(SphericalPoints);
	%% allocation
	Field=zeros(npoints,3,2);
	%% getting z-coordinates
	R=SphericalPoints(:,1);
	Theta=SphericalPoints(:,2);
	Phi=SphericalPoints(:,3);
	%% getting indices k and l
	k=Ind(1);
	l=Ind(2);
	%% function handle for associated legendre polynomials
	P=@(m,n,x) getAssoLegendre([m;n],x);
	%% cosine case
	%% e_r direction
	% empty
	%% e_theta direction
	if not(k==0) % avoid division by zero / else zero
		Field(:,2,1)=sin(k*Phi).*(P(k+1,l-1,cos(Theta)) ...
								+(k+l-1)*(k+l)*P(k-1,l-1,cos(Theta)))/2/k;
	end
	%% e_phi direction
	Field(:,3,1)=-cos(k*Phi).*(P(k+1,l,cos(Theta))...
								-(k+l)*(l-k+1)*P(k-1,l,cos(Theta)))/2;
	%% if radial function provided multiply by radial function
	if nargin==3
		Field(:,2,1)=f_r_cos(R).*Field(:,2,1);
		Field(:,3,1)=f_r_cos(R).*Field(:,3,1);
	end
	%% sin case
	%% e_r direction
	% empty
	%% e_theta direction
	if not(k==0) % avoid division by zero / else zero
		Field(:,2,2)=cos(k*Phi).*(P(k+1,l-1,cos(Theta)) ...
								+(k+l-1)*(k+l)*P(k-1,l-1,cos(Theta)))/2/k;
	end
	%% e_phi direction
	Field(:,3,2)=-sin(k*Phi).*(P(k+1,l,cos(Theta))...
								-(k+l)*(l-k+1)*P(k-1,l,cos(Theta)))/2;
	%% if radial function provided multiply by radial function
	if nargin==3
		Field(:,2,2)=f_r_sin(R).*Field(:,2,2);
		Field(:,3,2)=f_r_sin(R).*Field(:,3,2);
	end
end