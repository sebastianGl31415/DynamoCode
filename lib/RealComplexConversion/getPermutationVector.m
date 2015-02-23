function [p]=getPermutationVector(a)
% getPermutationVector -- returns a vector p in order to reorder row of 
% identity matrix I in order to get permutation matrix P.
%
% [p]=getPermutationVector(a)
%
% The permutation vector p is generated using a sequence starting from p(1)=a
% p(2)=2*a, p(3)=a+1, p(4)=a-1 iterating
%
%               / p(n-2)-1  if n > 3, n odd
%       p(n) = |                            .
%               \ p(n-2)+1  if n > 3, n odd 
%
	p=zeros(2*a,1);
	p(1)=a;
	p(2)=2*a;
	p(3:2:end-1)=[a+1:1:2*a-1]';
	p(4:2:end)=[a-1:-1:1]';
end