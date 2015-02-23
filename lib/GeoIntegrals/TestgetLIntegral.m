clear all;
close all;
compareMathematica=logical(1);
addpath('./../WignerSymbols/')
global params
global velocity
params.nmax=4;
velocity.k=[-2 0 2];
velocity.l=[0 1 2];
nmax=params.nmax;
L=zeros((nmax+1)^4,1);
Ind=zeros((nmax+1)^4,6);
lRow=0;
for n=0:nmax
	for m=-n:n
		%% getting all L-integral values for current pair
		[subL,subIndL]=getLIntegral(m,n);
		%% if there exists an index combination
		if not(isempty(subL))
			%% store integral values in L
			L(lRow+1:lRow+length(subL))=subL;
			%% store indices in IndK
			Ind(lRow+1:lRow+length(subL),:)=subIndL;
			%% increase row counter
			lRow=lRow+length(subL);
		end
	end
end
L=L(1:lRow);
Ind=Ind(1:lRow,:);
if compareMathematica==logical(1)
	load LIntegralMathematica.mat
	LMathematica=Expression1;
	clear Expression1;
	disp(['maximum error =  ',num2str(max(abs(L-LMathematica)))])
	figure
	hold on
	plot(imag(L),'k')
	plot(imag(LMathematica),'rx')
	xlabel('row of vector')
	ylabel('Integral L value')
	legend('octave/MATLAB-result','mathematica result')
else
	fname='InputList.mat';
	save('-mat-binary',fname,'Ind')
end