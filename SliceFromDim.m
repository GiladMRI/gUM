function Out=SliceFromDim(In,Idxs,Dim)
S.subs = repmat({':'},1,ndims(In));
S.subs{Dim} = Idxs;
S.type = '()';
Out=subsref(In,S);