function [O]=getO(k,l,Rh)
% getO --
	global velocity
	%% possible values for k and l
	K=[velocity.o{:,2}]';
	L=[velocity.o{:,3}]';
	%% find row corresponding to (k,l) mode
	[row,~,~]=find((K==k) & (L==l));
	if not(isempty(row))
		%% get function handle corresponding to (k,l) mode
		funcO=velocity.o{row,1};
		if (not(isvector(Rh)) && not(isscalar(Rh))) % matrix case
			%% get diagonal of sparse matrix
			rh=full(diag(Rh));
			n=length(rh);
			%% evaluate function
			O=funcO(rh);
			%% create sparse diagonal matrix for output
			O=sparse([1:n],[1:n],O,n,n);
		else % vector case
			O=funcO(Rh);
		end
	else
		O=0.0*Rh;
	end
end
