function [D]=applyRobinBoundary(D,options)
	%% loading global structure
	global params
	%% loading discretization variables
	diffMethod=params.diffMethod;
	h=params.h;
	%% input check for options
	if not(isstruct(options))
		error('2nd input must be a structure.')
	elseif not(any(strcmp('derivative',fieldnames(options))))
		error('derivative is not specified in 2nd input.')
	elseif not(any(strcmp('funcRobinBoundary',fieldnames(options))))
		error('funcRobinBoundary is not specified in 2nd input.')
	elseif not(any(strcmp('locationRobinBoundary',fieldnames(options))))
		error('locationRobinBoundary is not specified in 2nd input.')
	elseif not(any(strcmp('pointRobinBoundary',fieldnames(options))))
		error('pointRobinBoundary is not specified in 2nd input.')
	end
	%% reading options structure
	location=options.locationRobinBoundary;
	point=options.pointRobinBoundary;
	func=options.funcRobinBoundary;
	derivative=options.derivative;
	%% checking diffMethod
	if strcmp(diffMethod,'symmetric2ndOrder')
		%% checking location of neumann boundary
		if strcmp(location,'right') % x = 1
			%% checking of 1st or 2nd derivative
			switch derivative
				case '1' % 1st derivative y'
					a=func(point);
					D(end,:)=0.0; % set last row to zero
					D(end,end)=a(2)/a(1); % set remaining element
				case '2' % 2nd derivative y''
					a=func(point);
					D(end,:)=0.0; % set last row to zero
					D(end,end-1)=2/h^2; % set remaining elements
					D(end,end)=1/h^2*(2*h*a(2)/a(1)-2);
				otherwise
					error('derivative must either be a 1 or 2.')
			end
		else
			error(['Other than robin boundary at right are not implemented', ...
			' yet.'])
		end
	else
		error('Other than 2nd order symmetric differences are not implemented yet.')
	end
end