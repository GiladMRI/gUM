function Out=grepmat(In,Reps,Dims)
R=ones(1,30);
R(Dims)=Reps;
Out=repmat(In,R);