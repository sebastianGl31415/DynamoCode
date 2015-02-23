function [Wigner,jField]=WignerRecurrence(RecurrenceInput)
%%----------------------------------------------------------------------------%%
%% input check
%%----------------------------------------------------------------------------%%
	if (ne(size(RecurrenceInput,1),1) || ne(size(RecurrenceInput,2),5))
		error('wrong input size');
	end
	%% read input
	j2=RecurrenceInput(1);
	j3=RecurrenceInput(2);
	m1=RecurrenceInput(3);
	m2=RecurrenceInput(4);
	m3=RecurrenceInput(5);
%%----------------------------------------------------------------------------%%
%% preparation
%%----------------------------------------------------------------------------%%
	%% determine bounds for iteration
	jmin=max(abs(j2-j3),abs(m2+m3));
	jmax=j2+j3;
	jnum=jmax-jmin+1;
	jField=[jmin:jmax]';
	%% shortcuts for function definition
	a=[(j2-j3)^2 (j2+j3+1)^2 (m2+m3)^2];
	b=[(m2+m3)*(j2*(j2+1)-j3*(j3+1)) (m2-m3)];
	%% function definitions
	A=@(k) sqrt((k.^2-a(1)).*(a(2)-k.^2).*(k.^2-a(3)));
	B=@(k) (2*k+1).*(b(1)-b(2)*k.*(k+1));
	X=@(k) k.*A(k+1);
	Y=@(k) B(k);
	Z=@(k) (k+1).*A(k);
	NormW3j=@(phi,jField) phi/sqrt(sum((2*jField+1).*phi.^2));
	jIndex=@(k) k-jmin+1;
	%% defining flags for division by zero
	flag1=0;
	flag2=0;
	%% defining scale factor
	scaleFactor=1.0e3;
	%% allocation
	phip=zeros(jnum,1);
	phim=zeros(jnum,1);
	rs=zeros(jnum,1);
%%----------------------------------------------------------------------------%%
%% only one term present
%%----------------------------------------------------------------------------%%
	if jnum==1
		Wigner=1.0/sqrt(2*jmin+1);
		%% apply sign condition
		%% avoid division by zero
		if ne(sign(Wigner(end)),0)
			%% normal case
			Wigner=(-1)^(j2-j3+m2+m3)/sign(Wigner(end)).*Wigner;
		else
			%% case for Winger(end)=0 -> sign(Wigner(end))=0
			%% just assuming the convention sign(0)=1
			Wigner=(-1)^(j2-j3+m2+m3).*Wigner;
		end
		return;
	end
%%----------------------------------------------------------------------------%%
%% forward iteration for lower non-classical region from jmin to jm
%%----------------------------------------------------------------------------%%
	if (((m1==0) && (m2==0)) && (m3==0)) % all m's are zero
		phim(jIndex(jmin))=1.0;
		phim(jIndex(jmin+1))=0.0;
		jm=jmin+1;
	elseif Y(jmin)==0.0
		if X(jmin)==0.0 % case of division by zero, 2nd term is undefined
			flag1=1;
			jm=jmin;
		else % 2nd term is zero
			phim(jIndex(jmin))=1.0;
			phim(jIndex(jmin+1))=0.0;
			jm=jmin+1;
		end
	elseif X(jmin)*Y(jmin)>0.0 % 2nd term outside non-classical region
 		phim(jIndex(jmin))=1.0;
 		phim(jIndex(jmin+1))=-Y(jmin)/X(jmin);
		jm=jmin;
	else
		rs(jIndex(jmin))=-X(jmin)/Y(jmin);
		jm=jmax;
		for k=jmin+1:jmax-1
			nom=X(k);
			denom=Y(k)+Z(k)*rs(jIndex(k-1));
			if (((abs(nom)>abs(denom)) || (nom*denom>=0.0)) ... % condition for
						|| (denom==0.0))															% classical region
				jm=k-1;																						% or division by zero
				break;
			else % non-classical region
				rs(jIndex(k))=-nom/denom;
			end
		end
		phim(jIndex(jm))=1.0;
		for k=1:jm-jmin % calculate the products
			phim(jIndex(jm-k))=phim(jIndex(jm-k+1))*rs(jIndex(jm-k));
		end
		if jm==jmin % calculate at least two terms for recurrence
			phim(jIndex(jmin+1))=-Y(jmin)/X(jmin);
			jm=jmin+1;
		end
	end
	if jm==jmax % all terms are calculated
		%% apply normalization
		Wigner=NormW3j(phim,jField);
		%% apply sign condition
		if ne(sign(Wigner(end)),0)
			%% normal case
			Wigner=(-1)^(j2-j3+m2+m3)/sign(Wigner(end)).*Wigner;
		else
			%% case for Winger(end)=0 -> sign(Wigner(end))=0
			%% just assuming the convention sign(0)=1
			Wigner=(-1)^(j2-j3+m2+m3).*Wigner;
		end
		return % done
	end
%%----------------------------------------------------------------------------%%
%% backward iteration for upper non-classical region from jmax to jp
%%----------------------------------------------------------------------------%%
	if (((m1==0) && (m2==0)) && (m3==0)) % all m's are zero
		phip(jIndex(jmax))=1.0;
		phip(jIndex(jmax-1))=0.0;
		jp=jmax-1;
	elseif Y(jmax)==0.0
		if Z(jmax)==0.0 % case of division by zero, 2nd term is undefined
			flag2=1;
			jp=jmax;
		else
			phip(jIndex(jmax))=1.0;
			phip(jIndex(jmax-1))=-Y(jmax)/Z(jmax);
			jp=jmax-1;
		end
	elseif Y(jmax)*Z(jmax)>0.0 % 2nd term outside non-classical region
		phip(jIndex(jmax))=1.0;
		phip(jIndex(jmax-1))=-Y(jmax)/Z(jmax);
		jp=jmax-1;
	else
		rs(jIndex(jmax))=-Z(jmax)/Y(jmax);
		jp=jmin;
		for k=jmax-1:-1:jm
			nom=Z(k);
			denom=Y(k)+X(k)*rs(jIndex(k+1));
			if (((abs(nom)>abs(denom)) || (nom*denom>=0.0)) ...	% condition for
						|| (denom==0.0))															% classical region
				jp=k+1;																						% or division by zero
				break;
			else % non-classical region
				rs(jIndex(k))=-nom/denom;
			end
		end
		phip(jIndex(jp))=1.0;
		for k=1:jmax-jp % calculate the products
			phip(jIndex(jp+k))=phip(jIndex(jp+k-1))*rs(jIndex(jp+k));
		end
		if jp==jmax % calculate at least two terms for recurrence
			phip(jIndex(jmax-1))=-Y(jmax)/Z(jmax);
			jp=jmax-1;
		end
	end
	if jp==jmin % all terms are calculated
		%% apply normalization
		Wigner=NormW3j(phip,jField);
		%% apply sign condition
		if ne(sign(Wigner(end)),0)
			%% normal case
			Wigner=(-1)^(j2-j3+m2+m3)/sign(Wigner(end)).*Wigner;
		else
			%% case for Winger(end)=0 -> sign(Wigner(end))=0
			%% just assuming the convention sign(0)=1
			Wigner=(-1)^(j2-j3+m2+m3).*Wigner;
		end
		return % done
	end
%%----------------------------------------------------------------------------%%
%% 3-term recurrence iteration for classical region from jm+1 to jp-1
%%----------------------------------------------------------------------------%%
	if ((flag1==0) && (flag2==0))
		jmid=floor((jm+jp)/2);
		%% forward 3-term recurrence
		pt=phim;
		if (((ne(jmin,jmid-1) && ne(jmin,jmid)) && ne(jm,jmin)) && ... 
						ne(phim(jIndex(jm)),0.0)) 
			for k=jm:jmid-1
				phim(jIndex(k+1))=-(Z(k)*phim(jIndex(k-1))+Y(k)*phim(jIndex(k)))/X(k);
				if abs(phim(jIndex(k+1)))>1.0 % avoid overflows
						phim=phim/scaleFactor;
				end
				if ((abs(phim(jIndex(k+1))/phim(jIndex(k-1)))<1.0) && ...
							ne(phim(jIndex(k+1)),0.0))
					jmid=k+1;	% avoid upward iteration for decreasing values
					break;		% use downward iteration instead
				end
			end
			if ne(phim(jIndex(jmid-1)),0.0) && ...
						(abs(phim(jIndex(jmid))/phim(jIndex(jmid-1))) < 1.0e-6)
																								% avoid that the midpoint value
				jmid=jmid-1;														% is not zero or close to it
			end
		elseif (((jm==jmin) && ne(jmin,jmid)) && ne(phim(jIndex(jm)),0.0))
			for k=jm+1:jmid-1
				phim(jIndex(k+1))=-(Z(k)*phim(jIndex(k-1))+Y(k)*phim(jIndex(k)))/X(k);
				if abs(phim(jIndex(k+1)))>1.0 % avoid overflows
						phim=phim/scaleFactor;
				end
				if ((abs(phim(jIndex(k+1))/phim(jIndex(k-1)))<1.0) && ...
					ne(phim(jIndex(k+1)),0.0))
					jmid=k+1;	% avoid upward iteration for decreasing values
					break;		% use downward iteration instead
				end
			end
		elseif phim(jIndex(jm))==0.0	% avoid divison by zero
																	% iterate one term forward
				phim(jIndex(jm+1))=-(Z(jm)*phim(jIndex(jm-1))...
															+Y(jm)*phim(jIndex(jm)))/X(jm);
				jmid=jm+1;
		end
		%% backward 3-term recurrence
		for k=jp:-1:jmid+1
			phip(jIndex(k-1))=-(X(k)*phip(jIndex(k+1))+Y(k)*phip(jIndex(k)))/Z(k);
			if abs(phip(jIndex(k-1)))>1.0 % avoid overflows
				phip=phip/scaleFactor;
			end
		end
		if jmid==jmax
				phi=phim;
		elseif jmid==jmin
				phi=phip;
		else
			phi=[phim(jIndex(jmin):jIndex(jmid))*phip(jIndex(jmid))...
							/phim(jIndex(jmid)); phip(jIndex(jmid+1):jIndex(jmax))];
		end
	elseif (flag1==1) && (flag2==0)	% iteration in downward direction only
		for k=jp:-1:jmin+1
			phip(jIndex(k-1))=-(X(k)*phip(jIndex(k+1))+Y(k)*phip(jIndex(k)))/Z(k);
			if abs(phip(jIndex(k-1)))>1.0 % avoid overflows
				phip(jIndex(k-1):jIndex(jmax))=phip(jIndex(k-1):jIndex(jmax))...
																				/scaleFactor;
			end
		end
		phi=phip;
	elseif (flag1==0) && (flag2==1)	% iteration in upward direction only
		for k=jn:jp-1
			phim(jIndex(k+1))=-(Z(k)*phim(jIndex(k-1))+Y(k)*phim(jIndex(k)))/X(k);
			if abs(phi(jIndex(k+1)))>1.0 % avoid overflows
				phim(jIndex(jmin):jIndex(k+1))=phim(jindex(jmin):jindex(k+1))...
																					/scaleFactor;
			end
		end
		phi=phim;
	elseif (flag1==1) && (flag2==1)
		error('fatal error in Wigner3j, flag1 and flag2 are set.')
	end
	%% apply normalization
	Wigner=NormW3j(phi,jField);
	%% apply sign condition
	if ne(sign(Wigner(end)),0)
		%% normal case
		Wigner=(-1)^(j2-j3+m2+m3)/sign(Wigner(end)).*Wigner;
	else
		%% case for Winger(end)=0 -> sign(Wigner(end))=0
		%% just assuming the convention sign(0)=1
		Wigner=(-1)^(j2-j3+m2+m3).*Wigner;
	end
end