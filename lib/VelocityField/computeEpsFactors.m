function [oCosOut,oSinOut,pCosOut,pSinOut]=computeEpsFactors(oCosIn,...
																							oSinIn,pCosIn,pSinIn)
	global velocity
	if strcmp(velocity.referenceCase,'Gubbins2000')
		if not(isfield(velocity,'M')) || not(isfield(velocity,'D'))
			error('M and D must be specified in velocity structure.')
		else
			M=velocity.M;
			D=velocity.D;
		end
		[oCosIn,pCosIn]=getVelocityNormFactors(oCosIn,pCosIn);
		[oSinIn,pSinIn]=getVelocityNormFactors(oSinIn,pSinIn);
		%% getting alpha, beta, gamma, delta / Gubbins equation~(2.5)
		tau=zeros(1,1);
		row=0;
		for k=1:size(oCosIn,1)
			if not(isempty(oCosIn(k,4)))
				tau(row+1)=oCosIn{k,4};
				row=row+1;
			end
		end
		for k=1:size(oSinIn,1)
			if not(isempty(oSinIn{k,4}))
				tau(row+1)=oSinIn{k,4};
				row=row+1;
			end
		end
		varpi=zeros(3,1);
		row=0;
		for k=1:size(pCosIn,1)
			if not(isempty(pCosIn{k,4}))
				varpi(row+1)=pCosIn{k,4};
				row=row+1;
			end
		end
		for k=1:size(pSinIn,1)
			if not(isempty(pSinIn{k,4}))
				varpi(row+1)=pSinIn{k,4};
				row=row+1;
			end
		end
		%% only applicable to Gubbins2000 case
		alpha=tau(1);
		betta=varpi(1);
		ggamma=varpi(2);
		delta=varpi(3);
		%% seeting eps00 / Gubbins equation~(2.6)
		eps00=sign(D)*sqrt(D/alpha);
		%% seeting eps01 / Gubbins equation~(2.7)
		eps01=sign(M)*sqrt(M/betta);
		%% kinetic energy should be unity
		eps02=sqrt((alpha*eps00^2+betta*eps01^2)/(ggamma+delta));
		%% in text before equation~(2.5)
		eps03=eps02;
		%% rewriting tau and varpi
		tau(1)=eps00;
		varpi(1)=eps01;
		varpi(2)=eps02;
		varpi(3)=eps03;
		%% preparing output
		oCosOut=cell(size(oCosIn,1),3);
		oSinOut=cell(size(oSinIn,1),3);
		pCosOut=cell(size(pCosIn,1),3);
		pSinOut=cell(size(pSinIn,1),3);
		oCosOut(:,2:3)=oCosIn(:,2:3);
		oSinOut(:,2:3)=oSinIn(:,2:3);
		pCosOut(:,2:3)=pCosIn(:,2:3);
		pSinOut(:,2:3)=pSinIn(:,2:3);
		%% modifying all functions handles
		row=0;
		for k=1:size(oCosIn,1)
			if not(isempty(oCosIn(k,4)))
				strfunc=func2str(oCosIn{k,1});
				pos=strfind(strfunc,'@(');
				strfunc=strcat(strfunc(1:pos-1),strfunc(pos+4:end));
				strnormfunc=strcat('oCosOut{k,1}=@(r)',num2str(tau(row+1),16),...
															'*(',strfunc,');');
				eval(strnormfunc);
				row=row+1;
			else
				oCosOut{k,1}=oCosIn{k,1};
			end
		end
		for k=1:size(oSinIn,1)
			if not(isempty(oSinIn{k,4}))
				strfunc=func2str(oSinIn{k,1});
				pos=strfind(strfunc,'@(');
				strfunc=strcat(strfunc(1:pos-1),strfunc(pos+4:end));
				strnormfunc=strcat('oSinOut{k,1}=@(r)',num2str(tau(row+1),16),...
															'*(',strfunc,');');
				eval(strnormfunc);
				row=row+1;
			else
				oSinOut{k,1}=oSinIn{k,1};
			end
		end
		row=0;
		for k=1:size(pCosIn,1)
			if not(isempty(pCosIn{k,4}))
				strfunc=func2str(pCosIn{k,1});
				pos=strfind(strfunc,'@(');
				strfunc=strcat(strfunc(1:pos-1),strfunc(pos+4:end));
				strnormfunc=strcat('pCosOut{k,1}=@(r)',num2str(varpi(row+1),16),...
															'*(',strfunc,');');
				eval(strnormfunc);
				row=row+1;
			else
				pCosOut{k,1}=pCosIn{k,1};
			end
		end
		for k=1:size(pSinIn,1)
			if not(isempty(pSinIn{k,4}))
				strfunc=func2str(pSinIn{k,1});
				pos=strfind(strfunc,'@(');
				strfunc=strcat(strfunc(1:pos-1),strfunc(pos+4:end));
				strnormfunc=strcat('pSinOut{k,1}=@(r)',num2str(varpi(row+1),16),...
															'*(',strfunc,');');
				eval(strnormfunc);
				row=row+1;
			else
				pSinOut{k,1}=pSinIn{k,1};
			end
		end
	elseif strcmp(velocity.referenceCase,'ToroidalTestCase')
		[oCosIn,pCosIn]=getVelocityNormFactors(oCosIn,pCosIn);
		[oSinIn,pSinIn]=getVelocityNormFactors(oSinIn,pSinIn);
		tau=zeros(1,1);
		row=0;
		for k=1:size(oCosIn,1)
			if not(isempty(oCosIn(k,4)))
				tau(row+1)=oCosIn{k,4};
				row=row+1;
			end
		end
		for k=1:size(oSinIn,1)
			if not(isempty(oSinIn{k,4}))
				tau(row+1)=oSinIn{k,4};
				row=row+1;
			end
		end
		for k=1:size(pCosIn,1)
			if not(isempty(pCosIn{k,4}))
				varpi(row+1)=pCosIn{k,4};
				row=row+1;
			end
		end
		for k=1:size(pSinIn,1)
			if not(isempty(pCosIn{k,4}))
				tau(row+1)=pCosIn{k,4};
				row=row+1;
			end
		end
		%% preparing output
		oCosOut=cell(size(oCosIn,1),3);
		oSinOut=cell(size(oSinIn,1),3);
		pCosOut=cell(size(pCosIn,1),3);
		pSinOut=cell(size(pSinIn,1),3);
		oCosOut(:,2:3)=oCosIn(:,2:3);
		oSinOut(:,2:3)=oSinIn(:,2:3);
		pCosOut(:,2:3)=pCosIn(:,2:3);
		pSinOut(:,2:3)=pSinIn(:,2:3);
		tau=1/sqrt(tau);
		%% modifying all functions handles
		row=0;
		for k=1:size(oCosIn,1)
			if not(isempty(oCosIn(k,4)))
				strfunc=func2str(oCosIn{k,1});
				pos=strfind(strfunc,'@(');
				strfunc=strcat(strfunc(1:pos-1),strfunc(pos+4:end));
				strnormfunc=strcat('oCosOut{k,1}=@(r)',num2str(tau(row+1),16),...
															'*(',strfunc,');');
				eval(strnormfunc);
				row=row+1;
			else
				oCosOut{k,1}=oCosIn{k,1};
			end
		end
		for k=1:size(oSinIn,1)
			if not(isempty(oSinIn{k,4}))
				strfunc=func2str(oSinIn{k,1});
				pos=strfind(strfunc,'@(');
				strfunc=strcat(strfunc(1:pos-1),strfunc(pos+4:end));
				strnormfunc=strcat('oSinOut{k,1}=@(r)',num2str(tau(row+1),16),...
															'*(',strfunc,');');
				eval(strnormfunc);
				row=row+1;
			else
				oSinOut{k,1}=oSinIn{k,1};
			end
		end
		row=0;
		for k=1:size(pCosIn,1)
			if not(isempty(pCosIn{k,4}))
				strfunc=func2str(pCosIn{k,1});
				pos=strfind(strfunc,'@(');
				strfunc=strcat(strfunc(1:pos-1),strfunc(pos+4:end));
				strnormfunc=strcat('pCosOut{k,1}=@(r)',num2str(tau(row+1),16),...
															'*(',strfunc,');');
				eval(strnormfunc);
				row=row+1;
			else
				pCosOut{k,1}=pCosIn{k,1};
			end
		end
		for k=1:size(pSinIn,1)
			if not(isempty(pSinIn{k,4}))
				strfunc=func2str(pSinIn{k,1});
				pos=strfind(strfunc,'@(');
				strfunc=strcat(strfunc(1:pos-1),strfunc(pos+4:end));
				strnormfunc=strcat('pSinOut{k,1}=@(r)',num2str(tau(row+1),16),...
															'*(',strfunc,');');
				eval(strnormfunc);
				row=row+1;
			else
				pSinOut{k,1}=pSinIn{k,1};
			end
		end
	end
end