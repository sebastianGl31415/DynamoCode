clear all;
close all;
compareMathematica=logical(1);
addpath('./../WignerSymbols/')
global params
global velocity
params.nmax=4;
velocity.k=[-2 0 2];
velocity.l=[0 1 2 3];
nmax=params.nmax;
K=zeros((nmax+1)^4,1);
Ind=zeros((nmax+1)^4,6);
kRow=0;
for n=0:nmax
	for m=-n:n
		%% getting all K-integral values for current pair
		[subK,subIndK]=getKIntegral(m,n);
		%% if there exists an index combination
		if not(isempty(subK))
			%% store integral values in K
			K(kRow+1:kRow+length(subK))=subK;
			%% store indices in IndK
			Ind(kRow+1:kRow+length(subK),:)=subIndK;
			%% increase row counter
			kRow=kRow+length(subK);
		end
	end
end
K=K(1:kRow);
Ind=Ind(1:kRow,:);
if compareMathematica==logical(1)
	load KIntegralMathematica.mat
	KMathematica=Expression1;
	clear Expression1;
	disp(['maximum error =  ',num2str(max(abs(K-KMathematica)))])
	figure
	hold on
	plot(K,'k')
	plot(KMathematica,'rx')
	xlabel('row of vector')
	ylabel('Integral K value')
	legend('octave/MATLAB-result','mathematica result')
else
	fname='InputList.mat';
	save('-mat-binary',fname,'Ind')
end