function [Tor Pol]=computeShericalSeries(funcIn,Indmax,separationCase)
% computeShericalSeries -- return cell arrays Tor and Pol containing cosine and 
%             sin function handles of the series expansion of the function 
%             funcIn in poloidal and toroidal vector fields.
%
%    [Tor Pol]=computeShericalSeries(funcIn,Indmax,separationCase)
%
% author: Sebastian Glane
	addpath('./../SpecialFunctions/')
	%% input check for 1st input
	if not(iscell(funcIn)) % test if cell array
		error('1st input must be cell array.')
	elseif not(size(funcIn,1)==3)
		error('1st input must be at least a 3-by-1 cell')
	else % test if function handles
		for k=1:size(funcIn,1)
			if not(isa(funcIn{k},'function_handle'))
				error('1st input must be cell array of functions handles .')
			end
		end
	end
	if strcmp(varargin{1},'r_thetaphi')
		r_thetaphi=logical(1);
		if not(size(funcIn,2)==2)
			error('1st input does not match with separation indicator.')
		end
	elseif strcmp(separationCase,'rtheta_phi')
		error('This separation case is not implemented yet.')
	elseif strcmp(separationCase,'rphi_theta')
		error('This separation case is not implemented yet.')
	elseif strcmp(separationCase,'r_phi_theta')
		error('This separation case is not implemented yet.')
	else
		error('3rd input does not match the convenient declarations.')
	end
	%% input check for 2nd and 3rd input
	if (not(isscalar(Indmax))||not(rem(Indmax,2)==0||rem(Indmax+1,2)==0))||...
			Indmax<=0
		error('2nd input must be positive scalar integer.')
	elseif not(ischar(separationCase))
		error('3rd input must be a string.')
	end
	%% getting e_r, e_theta and e_phi function handles of vector field
	if r_thetaphi
		f_er=funcIn(1,1:2)';
		f_theta=funcIn(2,1:2)';
		f_phi=funcIn(3,1:2)';
	end
	rmpath('./../SpecialFunctions/')
end