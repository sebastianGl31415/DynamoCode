function [P]=getP(k,l,Rh)
% getP --
	global velocity
	%% possible values for k and l
	K=[velocity.p{:,2}]';
	L=[velocity.p{:,3}]';
	%% find row corresponding to (k,l) mode
	[row,~,~]=find((K==k) & (L==l));
	if not(isempty(row))
		%% get function handle corresponding to (k,l) mode
		funcP=velocity.p{row,1};
		if (not(isvector(Rh)) && not(isscalar(Rh))) % matrix case
			%% get diagonal of sparse matrix
			rh=full(diag(Rh));
			n=length(rh);
			%% evaluate function
			P=funcP(rh);
			%% create sparse diagonal matrix for output
			P=sparse([1:n],[1:n],P,n,n);
		else % vector case
			P=funcP(Rh);
		end
	else
		P=0.0*Rh;
	end
end
