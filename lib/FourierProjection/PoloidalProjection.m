function [PolProjection]=PoloidalProjection(funcIn,Ind,varargin)
	%% check if GNU Octave or MATLAB
	if exist('OCTAVE_VERSION')
		isOctave=logical(1);
	else
		isOctave=logical(0);
	end
	%% input check for Ind
	if not(size(Ind,1)==2) || not(size(Ind,2)==2)
		error('1st input must a vector of length 2.')
	elseif Ind(1)>Ind(2)
		error('1st index is larger than 2nd.')
	else 
		k=Ind(1);
		l=Ind(2);
	end
	%% default values for separation indicator
	rthetaphi=logical(1);
	r_thetaphi=logical(0);
	rtheta_phi=logical(0);
	rphi_theta=logical(0);
	r_phi_theta=logical(0);
	%% input check
	if not(iscell(funcIn)) % test if cell array
		error('1st input must be cell array.')
	else % test if function handles
		for k=1:size(funcIn)
			if not(isa(funcIn{k},'function_handle'))
				error('1st input must be cell array of functions handles .')
			end
		end
	end
	%% varargin input check
	if nargin==3
		if not(ischar(varargin{1})) % test if string
			error('If called with 3 inputs 3rd input must be string.')
		end
		rthetaphi=logical(0);
		%% checking for separation indicator
		if strcmp(varargin{1},'r_thetaphi')
			r_thetaphi=logical(1);
			if size(funcIn,1)==2 && size(funcIn,2)==2
				f_r=funcIn{:,1};
				f_thetaphi=funcIn{:,2};
			else
				error('1st input does not match with separation indicator.')
			end
		elseif strcmp(varargin{1},'rtheta_phi')
			rtheta_phi=logical(1);
			if size(funcIn,1)==2
				f_rtheta=funcIn{:,1};
				f_phi=funcIn{:,2};
			else
				error('1st input does not match with separation indicator.')
			end
		elseif strcmp(varargin{1},'rphi_theta')
			rphi_theta=logical(1);
			if size(funcIn,1)==2 && size(funcIn,2)==2
				f_rphi=funcIn{:,1};
				f_theta=funcIn{:,2};
			else
				error('1st input does not match with separation indicator.')
			end
		elseif strcmp(varargin{1},'r_phi_theta')
			r_phi_theta=logical(1);
			if size(funcIn,1)==3 && size(funcIn,2)==2
				f_r=funcIn{:,1};
				f_theta=funcIn{:,2};
				f_phi=funcIn{:,3};
			else
				error('1st input does not match with separation indicator.')
			end
			else
			error('3rd input does not match the convenient declarations.')
		end
	elseif nargin==4
		%% checking for discrete r-grid
		if not(ischar(varargin{1})) % test if string
			error('If called with 4 inputs 3rd input must be string.')
		end
		if not(isscalar(varargin{2})) && not(ismatrix(varargin{2})) % test if vector
				rh=varargin{2};
		else
			error('If called with 4 inputs 4th input must be a vector.')
		end
		rthetaphi=logical(0);
		%% checking for separation indicator
		if strcmp(varargin{1},'r_thetaphi')
			r_thetaphi=logical(1);
			if size(funcIn,1)==2 && size(funcIn,2)==2
				f_r=funcIn{:,1};
				f_thetaphi=funcIn{:,2};
			else
				error('1st input does not match with separation indicator.')
			end
		elseif strcmp(varargin{1},'rtheta_phi')
			rtheta_phi=logical(1);
			if size(funcIn,1)==2
				f_rtheta=funcIn{:,1};
				f_phi=funcIn{:,2};
			else
				error('1st input does not match with separation indicator.')
			end
		elseif strcmp(varargin{1},'rphi_theta')
			rphi_theta=logical(1);
			if size(funcIn,1)==2 && size(funcIn,2)==2
				f_rphi=funcIn{:,1};
				f_theta=funcIn{:,2};
			else
				error('1st input does not match with separation indicator.')
			end
		elseif strcmp(varargin{1},'r_phi_theta')
			r_phi_theta=logical(1);
			if size(funcIn,1)==3 && size(funcIn,2)==2
				f_r=funcIn{:,1};
				f_theta=funcIn{:,2};
				f_phi=funcIn{:,3};
			else
				error('1st input does not match with separation indicator.')
		end
		error('The case of a discrete r-grid is not implemented yet.')
	else
		f_rthetaphi=funcIn{1};
	end
	%% integrator tolerance
	TOL=1e-12;
	%% getting normalization constant
	N=getNormalizationConstant(k,l);
	if rthetaphi
		error('This separation case is not implemented yet.')
	elseif r_thetaphi % integration over theta and phi at once
		if not(k==0)
			%% convert radial function handle to string
			str_func_r=func2str(f_r);
			pos_f_r=strfind(str_func_r,'@(');
			str_func_r=strcat(str_func_r(1:pos_f_r-1),str_func_r(pos_f_r+4:end));
			%% functions for integration sin and cosine part
			cos_int_func=@(theta,phi) f_thetaphi(theta,phi).*,... 
												;
			sin_int_func=@(theta,phi) f_thetaphi(theta,phi).*,... 
												;
			%% integration sin and cosine part
			if isOctave % GNU Octave case
				cos_int=dblquad(cos_int_func,0,pi,0,2*pi,TOL);
				sin_int=dblquad(sin_int_func,0,pi,0,2*pi,TOL);
			else % MATLAB case
				cos_int=integral2(cos_int_func,0,pi,0,2*pi,'AbsTol',TOL);
				sin_int=integral2(sin_int_func,0,pi,0,2*pi,'AbsTol',TOL);
			end
			%% rescaling sin and cosine part
			a_factor=2*N^2/l/(l+1)*cos_int;
			b_factor=2*N^2/l/(l+1)*sin_int;
			%% create string to define the cosine function handle
			str_func_r_cos=strcat('a = @(r) a_factor * r .* ',str_func_r,';');
			eval(str_func_r_cos);
			%% create string to define the cosine function handle
			str_func_r_sin=strcat('b = @(r) b_factor * r .* ',str_func_r,';');
			eval(str_func_r_sin);
		elseif k==0
			%% convert radial function handle to string
			str_func_r=func2str(f_r);
			pos_f_r=strfind(str_func_r,'@(');
			str_func_r=strcat(str_func_r(1:pos_f_r-1),str_func_r(pos_f_r+4:end));
			%% functions for integration sin and cosine part
			cos_int_func=@(theta,phi) f_thetaphi(theta,phi).*,... 
												CosSphericalHarmonics([k; l],theta,phi);
			%% integration sin and cosine part
			if isOctave % GNU Octave case
				cos_int=dblquad(cos_int_func,0,pi,0,2*pi,TOL);
			else % MATLAB case
				cos_int=integral2(cos_int_func,0,pi,0,2*pi,'AbsTol',TOL);
			end
			%% rescaling sin and cosine part
			a_factor=2*N^2/l/(l+1)*cos_int;
			%% create string to define the cosine function handle
			str_func_r_cos=strcat('a = @(r) a_factor * r .* ',str_func_r,';');
			eval(str_func_r_cos);
			b=@(r) 0*r;
		end
	elseif rtheta_phi
		error('This separation case is not implemented yet.')
	elseif rphi_theta
		error('This separation case is not implemented yet.')
	elseif r_phi_theta
	end
	PolProjection=[a b];
end