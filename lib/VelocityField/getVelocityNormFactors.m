function [outTor,outPol]=getVelocityNormFactors(inTor,inPol)
	if not(iscell(inTor)) || not(iscell(inTor))
		error('1st or 2nd input is not a cell array.')
	elseif not(size(inTor,2)==7) || not(size(inPol,2)==7)
		error('cell array must be n-by-6 arrays.')
	end
	outTor=cell(size(inTor,1),4);
	outTor(:,1:3)=inTor(:,1:3);
	outPol=cell(size(inPol,1),4);
	outPol(:,1:3)=inPol(:,1:3);
	for row=1:size(inTor,1)
		if not(isempty(inTor{row,5}))
			k=inTor{row,2};
			l=inTor{row,3};
			N=getNormalizationConstant(k,l);
			if isempty(inTor{row,6})
				polyVec=inTor{row,5};
				%% integrate squared polynomial part multiply by r^2
				intpolyVec=polyint(conv(conv(polyVec,polyVec),[1 0 0]));
				%% evaluate integral from 0 to 1
				NormTor=l*(l+1)/N^2/2;
				outTor{row,4}=(polyval(intpolyVec,1)-polyval(intpolyVec,0))*NormTor;
				if k==0
					outTor{row,4}=2*outTor{row,4};
				end
			else
				%% % convert function handle to string
				strnonPoly=func2str(inTor{row,6});
				%% find the sub string @( in function handle strings
				posnonPoly=strfind(nonPoly,'@(');
				%% kick out the string @(?? out of the function handle strings
				strnonPoly=strcat(strnonPoly(1:posnonPoly-1),...
											strnonPoly(posnonPoly+4:end));
				%% polynomial part
				polyVec=inTor{row,5};
				%% square polynomial part and multiply by r^2
				polyVec2Int=conv(conv(polyVec,polyVec),[1 0 0]);
				%% polynomial vector to string
				strVec2Int=num2str(polyVec2Int);
				%% create string to define the function handles for integration
				strdef_intfunc=strcat('intfunc=@(r) polyval([',strVec2Int,'],r).*',...
					strnonPoly,'.*',strnonPoly,';');
				%% create function handle for integration
				eval(strdef_intfunc);
				%% quadrature integration
				NormTor=l*(l+1)/N^2/2;
				intVal=quad(intfunc,0,1);
				outTor{row,4}=intVal*NormTor;
				if k==0
					outTor{row,4}=2*outTor{row,4};
				end
			end
		end
	end
	for row=1:size(inPol,1)
		if not(isempty(inPol{row,5}))
			k=inPol{row,2};
			l=inPol{row,3};
			N=getNormalizationConstant(k,l);
			if isempty(inPol{row,6})
				polyVec=inPol{row,5};
				%% integrate squared polynomial part (e_r part)
				intpolyVec01=polyint(conv(polyVec,polyVec));
				%% test division of polynomial by r
				[~,remainder]=deconv(polyVec,[1 0]);
				if not(all(remainder==0))% check for singularities
					error('singularity in polynomial.')
				else
					%% apply DOperator to polynomial part
					polyVecDOperator=polyder(polyVec)+deconv(polyVec,[1 0]);
				end
				intpolyVec02=polyint(conv(conv(polyVecDOperator,polyVecDOperator),...
															[1 0 0]));
				%% evaluate integrals from 0 to 1
				intVal01=polyval(intpolyVec01,1)+polyval(intpolyVec01,0);
				intVal02=polyval(intpolyVec02,1)-polyval(intpolyVec02,0);
				NormPol=l*(l+1)/N^2/2;
				outPol{row,4}=NormPol*(intVal01*l*(l+1)+intVal02);
				if k==0
					outPol{row,4}=2*outPol{row,4};
				end
			else
				%% convert function handle to string
				nonPoly=func2str(inPol{row,6});
				%% find the sub string @( in function handle strings
				posnonPoly=strfind(nonPoly,'@(');
				%% kick out the string @(?? out of the function handle strings
				nonPoly=strcat(nonPoly(1:posnonPoly-1),nonPoly(posnonPoly+4:end));
				%% polynomial part
				polyVec=inPol{row,5};
				%% square polynomial part and multiply by r^2
				polyVec2Int01=conv(polyVec,polyVec);
				%% polynomial vector to string
				strVec2Int01=num2str(polyVec2Int01);
				%% create string to define the function handles for integration
				strdef_intfunc01=strcat('intfunc01=@(r) polyval([',strVec2Int01,...
					'],r).*',nonPoly,'.*',nonPoly,';');
				%% create function handle for integration
				eval(strdef_intfunc01)
				%% quadrature intergration
				intVal01=quad(intfunc01,0,1);
				%% convert function handle to string
				strderivnonPoly=func2str(inPol{row,7});
				%% find the sub string @( in function handle strings
				posnonPoly=strfind(strderivnonPoly,'@(');
				%% kick out the string @(?? out of the function handle strings
				strderivnonPoly=strcat(strderivnonPoly(1:posnonPoly-1),...
														strderivnonPoly(posnonPoly+4:end));
				%% derivative of polynomial part multiplied by r
				polyVec2Int02=conv(polyder(polyVec),[1 0]);
				strVec2Int02=num2str(polyVec2Int02);
				%% polynomial part multiplied by r
				strVec2Int03=num2str(conv(polyVec,[1 0]));
				%% polynomial part divided by r multiplied by r
				polyVec2Int04=polyVec;
				strVec2Int04=num2str(polyVec2Int04);
				%% create string to define the function handles for integration
				strdef_intfunc02=strcat('intfunc02=@(r) (',...
						'polyval([',strVec2Int02,'],r).*',nonPoly,...
						' + polyval([',strVec2Int03,'],r).*(',strderivnonPoly,')',...
						' + polyval([',strVec2Int04,'],r).*',nonPoly,').^2;');
				%% create function handle for integration
				eval(strdef_intfunc02);
				%% quadrature integration
				intVal02=quad(intfunc02,0,1);
				NormPol=l*(l+1)/N^2/2;
				outPol{row,4}=NormPol*(l*(l+1)*intVal01+intVal02);
				if k==0
					outPol{row,4}=2*outPol{row,4};
				end
			end
		end
	end
end