function Data=cleanDataset(Data,delRows)
% cleanDataset -- removes points doubling from dataset to be written.
%
%  [Data]=cleanDataset(Data,delRows)
%
% author: Sebastian Glane

	%% get dimensions
	dimensions=Data.dimensions;
	nr=dimensions(1);
	ntheta=dimensions(2);
	nphi=dimensions(3);
	%% remove rows from all scalars
	if isfield(Data,'Scalars')
		Scalars=Data.Scalars;
		for k=1:size(Scalars,1)
			if not(isempty(Scalars{k,1}))
				S=reshape(Scalars{k,1},[prod(dimensions)+nr*ntheta 1]);
				S(delRows)=[];
				Scalars{k}=S;
			end
		end
		Data.Scalars=Scalars;
	end
	%% remove rows from all vectors
	if isfield(Data,'Vectors')
		Vectors=Data.Vectors;
		for k=1:size(Vectors,1)
			if not(isempty(Vectors{k,1}))
				V=reshape(Vectors{k,1},[prod(dimensions)+nr*ntheta 3]);
				V(delRows,:)=[];
				Vectors{k}=V;
			end
		end
		Data.Vectors=Vectors;
	end
	%% remove rows from all tensors
	if isfield(Data,'Tensors')
		Tensors=Data.Tensors;
		for k=1:size(Tensors,1)
			if not(isempty(Tensors{k,1}))
				T=reshape(Tensors{k,1},[prod(dimensions)+nr*ntheta 3 3]);
				T(delRows,:,:)=[];
				Tensors{k}=T;
			end
		end
		Data.Tensors=Tensors;
	end
end
