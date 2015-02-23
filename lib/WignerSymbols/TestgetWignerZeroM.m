clear all;
close all;
compareMathematica=logical(1);
global params
global velocity
params.nmax=40;
velocity.k=[-2 1 0 1 2];
velocity.l=[0 1 2 3 4 5];
nmax=params.nmax;
WignerZeroM=zeros((nmax+1)^4,1);
WignerZeroMPlus=zeros((nmax+1)^4,1);
Ind=zeros((nmax+1)^4,6);
IndPlus=zeros((nmax+1)^4,6);
wRow=0;
wpRow=0;
for n=0:nmax
	%% getting all L-integral values for current pair
	[subWZM,subIndWZM]=getWignerZeroM(n);
	%% if there exists an index combination
	if not(isempty(subWZM))
		%% store integral values in L
		WignerZeroM(wRow+1:wRow+length(subWZM))=subWZM;
		%% store indices in IndK
		Ind(wRow+1:wRow+length(subWZM),:)=subIndWZM;
		%% increase row counter
		wRow=wRow+length(subWZM);
	end
	%% getting all L-integral values for current pair
	[subWZMP,subIndWZMP]=getWignerZeroM(n,'+');
	%% if there exists an index combination
	if not(isempty(subWZMP))
		%% store integral values in L
		WignerZeroMPlus(wpRow+1:wpRow+length(subWZMP))=subWZMP;
		%% store indices in IndK
		IndPlus(wpRow+1:wpRow+length(subWZMP),:)=subIndWZMP;
		%% increase row counter
		wpRow=wpRow+length(subWZMP);
	end
end
WignerZeroM=WignerZeroM(1:wRow);
WignerZeroMPlus=WignerZeroMPlus(1:wpRow);
Ind=Ind(1:wRow,:);
IndPlus=IndPlus(1:wpRow,:);
if compareMathematica==logical(1)
	load WignerMathematica.mat
	WignerMathematica=Expression1;
	clear Expression1;
	load WignerPlusMathematica.mat
	WignerPlusMathematica=Expression1;
	clear Expression1;
	disp(['maximum error(standard)=',...
				num2str(max(abs(WignerZeroM-WignerMathematica)))])
	disp(['maximum error (plus)=',...
				num2str(max(abs(WignerZeroMPlus-WignerPlusMathematica)))])
	figure
	hold on
	plot(WignerZeroM,'k')
	plot(WignerMathematica,'rx')
	xlabel('row of vector')
	ylabel('ThreeJ value (standard)')
	legend('octave/MATLAB-result','mathematica result')
	figure
	hold on
	plot(WignerZeroMPlus,'k')
	plot(WignerPlusMathematica,'rx')
	xlabel('row of vector')
	ylabel('ThreeJ value (plus)')
	legend('octave/MATLAB-result','mathematica result')
else
	fname='InputList.mat';
	save('-mat-binary',fname,'Ind')
	fname='InputListPlus.mat';
	save('-mat-binary',fname,'IndPlus')
end
