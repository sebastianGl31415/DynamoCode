function []=displayOutput()
	if exist('OCTAVE_VERSION')
		fflush(stdout);
	else
		drawnow('update')
	end
end