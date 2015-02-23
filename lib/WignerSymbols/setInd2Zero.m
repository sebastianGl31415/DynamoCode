function [varargout]=setInd2Zero(Ind,varargin)
	for k=1:nargin
		varargin{k}(Ind)=0;
	end
	varargout=varargin;
end
