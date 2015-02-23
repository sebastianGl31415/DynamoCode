function []=initializeVelocity()
	global velocity;
	if not(isfield(velocity,'referenceCase'))
		error('referenceCase is not specified in velocity structure.')
	elseif not(ischar(velocity.referenceCase))
		error('referenceCasein velocity structure is not a character array.')
	end
	fprintf(['------- velocity field information -------------------------\n',...
'  velocity reference case: %s \n'],velocity.referenceCase);
	displayOutput()
	%% checking for velocity case
	if strcmp(velocity.referenceCase,'Gubbins2000')
		%% checking if p exists
		if isfield(velocity,'p') && isfield(velocity,'M') && isfield(velocity,'D')
			pp=velocity.p;
			fprintf(['  with the following parameter set:\n',...
'    -- p = %3.1f \n'],velocity.p);
			fprintf('    -- M = %4.3f\n',velocity.M);
			fprintf('    -- D = %4.3f\n',velocity.D);
		else
			error('For velocity case Gubbins2000 p, M and D needs to be specified.')
		end
		displayOutput()
		%% real valued toroidal fields including phantom sin modes
		oSin=cell(1,7);
		oCos=cell(1,7);
		%% k and l values for toroidal velocity modes
		oSin(1,2:3)=num2cell([0 1]); % phantom mode
		oCos(1,2:3)=num2cell([0 1]);
		%% function handles for toroidal velocity modes
		oSin{:,1}= @(r) 0*r;	% phantom mode
		oCos{:,1}= @(r) r.^2.*(1-r.^2);
		%% seperation of function
		%% in polynomial part and
		oSin{:,5}= [];	% phantom mode
		oCos{:,5}= conv([1 0 0],[-1 0 1]);
		%% real valued poloidal fields including phantom sin modes
		pSin=cell(2,7);
		pCos=cell(2,7);
		%% k and l values for poloidal velocity modes
		pSin(:,2)=num2cell([0;2]);
		pSin(:,3)=num2cell([2;2]);
		pCos(:,2)=num2cell([0;2]);
		pCos(:,3)=num2cell([2;2]);
		%% work around to get function hand with pp as a number
		strfunc01=['@(r) r.^4.*(1-r.^2).^2.*sin(',num2str(pp,16),'*pi*r)'];
		strfunc02=['@(r) r.^4.*(1-r.^2).^2.*cos(',num2str(pp,16),'*pi*r)'];
		%% function handles for poloidal velocity modes
		pSin(:,1)={@(r) 0*r; eval(strfunc01)};
		pCos(:,1)={@(r) r.^6.*(1-r.^2).^3; eval(strfunc02)};
		%% seperation of function
		%% in polynomial part
		pSin(:,5)={[];...	% phantom mode
								conv([1 0 0 0 0],conv([-1 0 1],[-1 0 1]))};
		pCos(:,5)={conv([1 0 0 0 0 0 0],conv(conv([-1 0 1],[-1 0 1]),[-1 0 1]));...
								conv([1 0 0 0 0],conv([-1 0 1],[-1 0 1]))};
		%% and non-polynomial part due to division by r including derivative
		strfunc01=['@(r) sin(',num2str(pp,16),'*pi*r)'];
		strfunc02=['@(r) cos(',num2str(pp,16),'*pi*r)'];
		pSin(:,6)={ []; eval(strfunc01)};
		pCos(:,6)={ []; eval(strfunc02)};
		strfunc01=['@(r) ',num2str(pp,16),'*pi*cos(',num2str(pp,16),'*pi*r)'];
		strfunc02=['@(r) -',num2str(pp,16),'*pi*sin(',num2str(pp,16),'*pi*r)'];
		pSin(:,7)={ []; eval(strfunc01)};
		pCos(:,7)={ []; eval(strfunc02)};
	elseif strcmp(velocity.referenceCase,'ToroidalTestCase')
		%% real valued poloidal fields including phantom sin modes
		oSin=cell(1,7);
		oCos=cell(1,7);
		%% k and l values for toroidal velocity modes
		oSin(1,2:3)=num2cell([0 1]); % phantom mode
		oCos(1,2:3)=num2cell([0 1]);
		%% function handles for toroidal velocity modes
		oSin{:,1}= @(r) 0*r;	% phantom mode
		oCos{:,1}= @(r) r.^2.*(1-r.^2);
		%% function handles for toroidal velocity modes
		oSin{:,1}= @(r) 0*r;	% phantom mode
		oCos{:,1}= @(r) r.^2.*(1-r.^2);
		%% seperation of function
		%% in polynomial part and
		oSin{:,5}= [];	% phantom mode
		oCos{:,5}= conv([1 0 0],[-1 0 1]);
		%% real valued poloidal fields including phantom sin modes
		pSin=cell(1,7);
		pCos=cell(1,7);
	else
		error('Other velocity cases are not implemented yet.')
	end
	%% normalization
	[oCos,oSin,pCos,pSin]=computeEpsFactors(oCos,oSin,pCos,pSin);
	%% converting to array from real to complex formulation
	[p]=convertReal2Complex(pSin,pCos);
	[o]=convertReal2Complex(oSin,oCos);
	%% writing to global velocity structure 
	velocity.p=p;
	velocity.o=o;
	%% getting k and l values as vectors
	k=[cell2mat(o(:,2)); cell2mat(p(:,2))];
	l=[cell2mat(o(:,3)); cell2mat(p(:,3))];
	%% no doubling of entries
	k=unique(k,'rows');
	l=unique(l,'rows');
	%% writing unique k and l values to global velocity structure
	velocity.k=k;	 % required to compute the K and L integrals
	velocity.l=l;
end