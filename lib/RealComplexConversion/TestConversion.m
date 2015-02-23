clear all
close all
clc
addpath('./../SpecialFunctions/');
global params
params.nmax=2;
params.nr=1;
nmax=params.nmax;
nr=params.nr;
%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
fprintf('-------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Cycling from RealModes to RealModesRev\n');
nmodes=(nmax+1)*(nmax+2);
RealModes=zeros(nmodes*nr,1);
IndRealModes=zeros(nmodes,2);
for l=0:nmax
	for k=0:l
		row=l*(l+1)+k*2+1;
		RealModes((row-1)*nr+1:row*nr)=randn(nr,1)+i*randn(nr,1);
		IndRealModes(row,1:2)=[k l];
	end
end
for l=1:nmax
	for k=1:l
		row=l*(l+1)+(k+1)*2;
		RealModes((row-1)*nr+1:row*nr)=randn(nr,1)+i*randn(nr,1);
		IndRealModes(row,1:2)=[k l];
	end
end
row=find(IndRealModes(3:end,2)==0)+2;
IndRealModes(row,:)=IndRealModes(row-1,:);
clear k  l nmodes row
%------------------------------------------------------------------------------%
[SinModes,CosModes]=getSinCosModes(RealModes,IndRealModes);
[RealModesRev,IndRealModesRev]=getRealModes(SinModes,CosModes);
%------------------------------------------------------------------------------%
fprintf('Comparing RealModes\n');
fprintf('any parts matching: %i\n',any(RealModes==RealModesRev));
fprintf('any real parts matching: %i\n', ... 
						any(real(RealModes)==real(RealModesRev)));
fprintf('any imaginary parts matching: %i\n', ... 
						any(imag(RealModes)==imag(RealModesRev)));
fprintf('all parts matching: %i\n',all(RealModes==RealModesRev));
fprintf('all real parts matching: %i\n', ... 
						all(real(RealModes)==real(RealModesRev)));
fprintf('all imaginary parts matching: %i\n', ... 
						all(imag(RealModes)==imag(RealModesRev)));
clear SinModes CosModes RealModes IndRealModes RealModesRev IndRealModesRev
%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
fprintf('-------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Cycling from SinModes/CosModes to SinModes/CosModes\n');
nsin=0.5*(nmax+1)*nmax;
ncos=0.5*(nmax+1)*(nmax+2);
SinModes=cell(nsin,3);
CosModes=cell(ncos,3);
rowcount=1;
for l=0:nmax
	for k=0:l
		CosModes(rowcount,1)=randn(nr,1)+i*randn(nr,1);
		CosModes(rowcount,2)=k;
		CosModes(rowcount,3)=l;
		rowcount=rowcount+1;
	end
end
rowcount=1;
for l=1:nmax
	for k=1:l
		SinModes(rowcount,1)=randn(nr,1)+i*randn(nr,1);
		SinModes(rowcount,2)=k;
		SinModes(rowcount,3)=l;
		rowcount=rowcount+1;
	end
end
clear k l nsin ncos rowcount
%------------------------------------------------------------------------------%
[RealModes,IndRealModes]=getRealModes(SinModes,CosModes);
[SinModesRev,CosModesRev]=getSinCosModes(RealModes,IndRealModes);
%------------------------------------------------------------------------------%
SinVal=cell2mat(SinModes(:,1));
SinValRev=cell2mat(SinModesRev(:,1));
CosVal=cell2mat(CosModes(:,1));
CosValRev=cell2mat(CosModesRev(:,1));
%------------------------------------------------------------------------------%
fprintf('Comparing SinModes\n');
fprintf('any parts matching: %i\n',any(SinVal==SinValRev));
fprintf('any real parts matching: %i\n', ... 
						any(real(SinVal)==real(SinValRev)));
fprintf('any imaginary parts matching: %i\n', ... 
						any(imag(SinVal)==imag(SinValRev)));
fprintf('all parts matching: %i\n',all(SinVal==SinValRev));
fprintf('all real parts matching: %i\n', ... 
						all(real(SinVal)==real(SinValRev)));
fprintf('all imaginary parts matching: %i\n', ... 
						all(imag(SinVal)==imag(SinValRev)));
%------------------------------------------------------------------------------%
fprintf('Comparing CosModes\n');
fprintf('any parts matching: %i\n',any(CosVal==CosValRev));
fprintf('any real parts matching: %i\n', ... 
						any(real(CosVal)==real(CosValRev)));
fprintf('any imaginary parts matching: %i\n', ... 
						any(imag(CosVal)==imag(CosValRev)));
fprintf('all parts matching: %i\n',all(CosVal==CosValRev));
fprintf('all real parts matching: %i\n', ... 
						all(real(CosVal)==real(CosValRev)));
fprintf('all imaginary parts matching: %i\n', ... 
						all(imag(CosVal)==imag(CosValRev)));
clear CosVal CosValRev SinVal SinValRev CosModes CosModesRev SinModes ... 
				SinModesRev RealModes IndRealModes
%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
fprintf('-------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Converting from ComplexModes to RealModes to ComplexModesRev\n');
%------------------------------------------------------------------------------%
ComplexModes=zeros((nmax+1)^2*nr,1);
for l=0:nmax
	for k=-l:-1
		row=l*(l+1)+k+1;
		ComplexModes((row-1)*nr+1:row*nr)=randn(nr,1)+i*randn(nr,1);
	end
	row=l*(l+1)+1;
	ComplexModes((row-1)*nr+1:row*nr)=randn(nr,1);
	for k=1:l
		roww=l*(l+1)+k+1;
		rowr=l*(l+1)-k+1;
		ComplexModes((roww-1)*nr+1:roww*nr)=conj(...
			ComplexModes((rowr-1)*nr+1:rowr*nr));
	end
end
clear k l row roww rowr
%------------------------------------------------------------------------------%
[RealModes,IndRealModes]=convertComplex2Real(ComplexModes);
[ComplexModesRev]=convertReal2Complex(RealModes,IndRealModes);
%------------------------------------------------------------------------------%
fprintf('Comparing ComplexModes\n');
fprintf('any parts matching: %i\n',... 
						any(abs(ComplexModesRev-ComplexModes)<10*eps));
fprintf('any real parts matching: %i\n', ... 
						any(real(ComplexModes-ComplexModesRev)<10*eps));
fprintf('any imaginary parts matching: %i\n', ... 
						any(imag(ComplexModes-ComplexModesRev)<10*eps));
fprintf('all parts matching: %i\n',...
						all(abs(ComplexModesRev-ComplexModes)<10*eps));
fprintf('all real parts matching: %i\n', ... 
						all(real(ComplexModes-ComplexModesRev)<10*eps));
fprintf('all imaginary parts matching: %i\n', ... 
						all(real(ComplexModes-ComplexModesRev)<10*eps));
clear ComplexModes ComplexModesRev RealModes IndRealModes
%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
fprintf('-------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Converting from RealModes to ComplexModes to to RealModesRev\n');
nmodes=(nmax+1)*(nmax+2);
RealModes=zeros(nmodes*nr,1);
IndRealModes=zeros(nmodes,2);
for l=0:nmax
	for k=0:l
		row=l*(l+1)+k*2+1;
		RealModes((row-1)*nr+1:row*nr)=randn(nr,1)+i*randn(nr,1);
		IndRealModes(row,1:2)=[k l];
	end
end
for l=1:nmax
	for k=1:l
		row=l*(l+1)+(k+1)*2;
		RealModes((row-1)*nr+1:row*nr)=randn(nr,1)+i*randn(nr,1);
		IndRealModes(row,1:2)=[k l];
	end
end
row=find(IndRealModes(3:end,2)==0)+2;
IndRealModes(row,:)=IndRealModes(row-1,:);
clear k  l nmodes row
%------------------------------------------------------------------------------%
[ComplexModes]=convertReal2Complex(RealModes,IndRealModes);
[RealModesRev,IndRealModesRev]=convertComplex2Real(ComplexModes);
%------------------------------------------------------------------------------%
fprintf('Comparing ComplexModes\n');
fprintf('any parts matching: %i\n',... 
						any(abs(RealModes-RealModesRev)<10*eps));
fprintf('any real parts matching: %i\n', ... 
						any(real(RealModes-RealModesRev)<10*eps));
fprintf('any imaginary parts matching: %i\n', ... 
						any(imag(RealModes-RealModesRev)<10*eps));
fprintf('all parts matching: %i\n',...
						all(abs(RealModes-RealModesRev)<10*eps));
fprintf('all real parts matching: %i\n', ... 
						all(real(RealModes-RealModesRev)<10*eps));
fprintf('all imaginary parts matching: %i\n', ... 
						all(real(RealModes-RealModesRev)<10*eps));
clear ComplexModes RealModesRev RealModes IndRealModes
%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
fprintf('-------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------\n');
rmpath('./../SpecialFunctions/');