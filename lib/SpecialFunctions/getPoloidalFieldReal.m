function [Field]=getPoloidalFieldReal(Ind,SphericalPoints,varargin)
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
			f_r_cos=funcIn{1};
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
	Field(:,1,1)=l*(l+1)*getCosSphericalHarmonic([k;l],Theta,Phi);
	%% e_theta direction
	Field(:,2,1)=cos(k*Phi).*((k+l)*(l-k+1)*P(k-1,l,cos(Theta)) ...
								-P(k+1,l,cos(Theta)))/2;
	%% e_phi direction
	if not(k==0) % avoid division by zero / else zero
		Field(:,3,1)=sin(k*Phi).*(P(k+1,l-1,cos(Theta)) ...
								+(k+l-1)*(k+l)*P(k-1,l-1,cos(Theta)))/2/k;
	end
	%% if radial function provided multiply by radial function
	if nargin==3
		%% e_r direction
		Field(:,2,1)=f_r_cos(R)./R.*Field(:,2,1);
		%% set all NaN and Inf to zero
		ind=find(bsxfun(@or,isnan(Field),isinf(Field)));
		Find(Ind)=0.0;
		%% DOperator / use polynomial routines deconv, polyderiv, polyval
		df_r=gradient(f_r_cos(R),R);
		f_over_r=f_r_cos(R)./R;
		%% set all NaN and Inf to zero
		ind=find(bsxfun(@or,isnan(f_over_r),isinf(f_over_r)));
		f_over_r(ind)=0.0;
		dOperator=df_r+f_over_r;
		%% e_theta / e_phi direction
		Field(:,2,1)=dOperator.*Field(:,2,1);
		Field(:,3,1)=dOperator.*Field(:,3,1);
	end
	%% sin case
	%% e_r direction
	Field(:,1,1)=l*(l+1)*getSinSphericalHarmonic([k;l],Theta,Phi);
	%% e_theta direction
	Field(:,2,1)=sin(k*Phi).*((k+l)*(l-k+1)*P(k-1,l,cos(Theta)) ...
								-P(k+1,l,cos(Theta)))/2;
	%% e_phi direction
	if not(k==0) % avoid division by zero / else zero
		Field(:,3,1)=cos(k*Phi).*(P(k+1,l-1,cos(Theta)) ...
								+(k+l-1)*(k+l)*P(k-1,l-1,cos(Theta)))/2/k;
	end
	%% if radial function provided multiply by radial function
	if nargin==3
		%% e_r direction
		Field(:,2,1)=f_r_sin(R)./R.*Field(:,2,1);
		%% set all NaN and Inf to zero
		ind=find(bsxfun(@or,isnan(Field),isinf(Field)));
		Find(Ind)=0.0;
		%% DOperator / use polynomial routines deconv, polyderiv, polyval
		df_r=gradient(f_r_sin(R),R);
		f_over_r=f_r_sin(R)./R;
		%% set all NaN and Inf to zero
		ind=find(bsxfun(@or,isnan(f_over_r),isinf(f_over_r)));
		f_over_r(ind)=0.0;
		dOperator=df_r+f_over_r;
		%% e_theta / e_phi direction
		Field(:,2,1)=dOperator.*Field(:,2,1);
		Field(:,3,1)=dOperator.*Field(:,3,1);
	end
end