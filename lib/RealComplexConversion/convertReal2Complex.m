function [Out]=convertReal2Complex(varargin)
% convertReal2Complex --
	%% check if MATLAB or GNU Octave
	if exist('OCTAVE_VERSION')
	else
		iscomplex=@(in) not(isreal(in));
	end
	%% input check
	if not(nargin==2)
		error('Function accepts exactly two input arguments.')
	elseif not(iscell(varargin{1}) && iscell(varargin{2})) &&...
			not(isvector(varargin{1}) && ismatrix(varargin{2}))
		error('Inputs must either be two cell arrays or a vector and a matrix.')
	end
	%% check which case to execute
	if iscell(varargin{1}) && iscell(varargin{2})
		%% convert input
		sinCell=varargin{1};
		cosCell=varargin{2};
		if not(all(size(sinCell)==size(cosCell)))
			error('Cell array do not equal in size.')
		end
		%% find rows if cos modes i.e. k = 0
		rowzero=find(cell2mat(sinCell(:,2))==0);
		%% find rows if mixed modes
		rownonzero=find(not(cell2mat(sinCell(:,2))==0));
		%% calculate number of modes
		nmodes=length(rowzero)+2*length(rownonzero);
		%% allocation of complex cell array
		compCell=cell(nmodes,3);
		%% check imaginary unit not be overwritten
		if not(iscomplex(i))
			error('i is not the imaginary unit anymore')
		end
		%% loop over mixed modes only
		for ind=1:length(rownonzero)
			%% retrieve infomation from real cell arrays
			funcSin=func2str(sinCell{rownonzero(ind),1});	% convert function handle to
			funcCos=func2str(cosCell{rownonzero(ind),1});	% string of the form
			k=sinCell{rownonzero(ind),2};									% @(x) fcos(x)/fsin(x)
			l=sinCell{rownonzero(ind),3};
			%% find the sub string @( in function handle strings
			posSin=strfind(funcSin,'@(');
			posCos=strfind(funcCos,'@(');
			%% kick out the string @(?? out of the function handle strings
			funcSin=strcat(funcSin(1:posSin-1),funcSin(posSin+4:end));
			funcCos=strcat(funcCos(1:posCos-1),funcCos(posCos+4:end));
			%% get normalization constant for real-complex conversion
			N=getNormalizationConstant(k,l);
			%% write indices in complex cell array
			compCell{(2*ind-1),2}=-k;	% negative k case
			compCell{2*ind,2}=k;			% positive k case
			compCell([(2*ind-1) 2*ind],3)={l; l};
			%% create string to define the complex function handles
			strfuncCompMinus=strcat('funcCompMinus = @(r) (',funcCos,...
				' + i * ',funcSin,' )/(2*N);');	% negative k case 
			strfuncCompPlus=strcat('funcCompPlus = @(r) (',funcCos,...
				' - i * ',funcSin,' )/(2*N);');	% positive k case 
			%% create complex function handles by executing the definition by string 
			eval(strfuncCompMinus);	% negative k case
			eval(strfuncCompPlus);	% positive k case
			%% writing function handles in cell array
			compCell{(2*ind-1),1}=funcCompMinus;
			compCell{2*ind,1}=funcCompPlus;
		end
		%% save rows written in complex cell array as shift for next loop
		rowshift=2*length(rownonzero);
		%% loop over cos modes only
		for ind=1:length(rowzero)
			%% retrieve infomation from real cell arrays
			funcCos=func2str(cosCell{rowzero(ind),1});	% convert function handle to
			k=cosCell{rowzero(ind),2};										% string of the form
			l=cosCell{rowzero(ind),3};										% @(x) fcos(x)/fsin(x)
			%% find the sub string @( in function handle strings
			posCos=strfind(funcCos,'@(');
			%% kick out the string @(?? out of the function handle strings
			funcCos=strcat(funcCos(1:posCos-1),funcCos(posCos+4:end));
			%% get normalization constant for real-complex conversion
			N=getNormalizationConstant(k,l);
			%% write indices in complex cell array
			compCell{rowshift+ind,2}=k;
			compCell{rowshift+ind,3}=l;
			%% create string to define the complex function handles
			strfuncComp=strcat('funcComp = @(r) (',funcCos,' )/N;');
			%% create complex function handles by executing the definition by string 
			eval(strfuncComp);
			%% writing function handles in cell array
			compCell{rowshift+ind,1}=funcComp;
		end
		Out=compCell;
	elseif (isvector(varargin{1}) && (not(isvector(varargin{2})) && ...
																					ismatrix(varargin{2})))
		%% convert modes
		RealModes=varargin{1};
		IndRealModes=varargin{2};
		global params
		nmax=params.nmax;
		nr=params.nr;
		%% number of real modes (observe that nmodes is always even)
		nmodes=(nmax+1)*(nmax+2);
		%% input check
		if not(size(RealModes,1)==nmodes*nr || ...
				size(RealModes,1)==nmodes*(nr+1))
			error('Length of the 1st input is not correct.')
		elseif not(size(RealModes,2)==1)
			error('Only vector inputs are accepted.')
		elseif not(size(IndRealModes,1)==nmodes)
			error('Length of the 2nd input is not correct.')
		end
		%% read k- and l-values
		k=IndRealModes(:,1);
		l=IndRealModes(:,2);
		%% normalization constant
		N=getNormalizationConstant(k,l);
		%% sparse diagonal matrix with every 2nd entry N^k_l on diagonal
		%% every 2nd entry since N^k_l = N^(-k)_l
		spdiagN=sparse([1:nmodes/2],[1:nmodes/2],1./(2*N(1:2:end)),... 
										nmodes/2,nmodes/2) ;
		%% check if poloidal or toroidal field by number of grid points
		if (size(RealModes,1)/nmodes)==nr
			np=nr;
		elseif (size(RealModes,1)/nmodes)==(nr+1)
			np=nr+1;
		else
			error('Something is wrong with the number of points.')
		end
		%% identity for grid points
		Igrid=speye(np);
		%% check imaginary unit not be overwritten
		if not(iscomplex(i))
			error('i is not the imaginary unit anymore')
		end
		%% matrix converting real pairs to complex pairs
		Mr2c=[1 -i;1 i];
		%% build large conversion matrix by kronecker products
		ConversionMatrix=kron(kron(spdiagN,Mr2c),Igrid);
		%% conversion by sparse matrix vector product
		ComplexModes=ConversionMatrix*RealModes;
		%% eliminate phantom modes
		[ComplexModes,Ind]=sortModeVector(ComplexModes,'red');
		row=find(Ind(:,1)==0);
		for k=1:length(row)
			r=row(k);
			ComplexModes((r-1)*np+1:r*np)=2*ComplexModes((r-1)*np+1:r*np);
		end
		Out=ComplexModes;
	else
		error('Something is wrong with the inputs.')
	end
end