clear all;
close all;
compareMathematica=logical(0);
global params
global velocity
params.nmax=40;
velocity.k=[-4 -3 -2 -1 0 -1 2 3 4];
velocity.l=[0 1 2 3 5 6];
nmax=params.nmax;
Wigner=zeros((nmax+1)^4,1);
Ind=zeros((nmax+1)^4,6);
wRow=0;
for n=0:nmax
	for m=-n:n
		%% getting all L-integral values for current pair
		[subWigner,subIndWigner]=getAllWigner(n,m);
		%% if there exists an index combination
		if not(isempty(subWigner))
			%% store integral values in L
			Wigner(wRow+1:wRow+length(subWigner))=subWigner;
			%% store indices in IndK
			Ind(wRow+1:wRow+length(subWigner),:)=subIndWigner;
			%% increase row counter
			wRow=wRow+length(subWigner);
		end
	end
end
Wigner=Wigner(1:wRow);
Ind=Ind(1:wRow,:);
if compareMathematica==logical(1)
	load WignerMathematica.mat
	WignerMathematica=Expression1;
	clear Expression1;
	disp(['maximum error =  ',num2str(max(abs(Wigner-WignerMathematica)))])
	figure
	hold on
	plot(Wigner,'k')
	plot(WignerMathematica,'rx')
	xlabel('row of vector')
	ylabel('ThreeJ value')
	legend('octave/MATLAB-result','mathematica result')
else
	fname='InputList.mat';
	save('-mat-binary',fname,'Ind')
end