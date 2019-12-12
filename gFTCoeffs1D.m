function Out=gFTCoeffs1D(Sz,Idxs)

if(mod(Sz,2)==1 && Sz>1)
%     Idxs=Idxs+0.5;
end

% Line= -(0:(Sz-1))*2*pi/Sz;
% Line=-linspaceWithHalfStep(-Sz/2,Sz/2,Sz)*2*pi/Sz;
Line=-linspaceWithHalfStep(-0.5,0.5,Sz)*2*pi;
% Line=-linspace(-Sz/2,Sz/2,Sz)*2*pi/Sz + 4*pi;
M=(Idxs(:))*Line;

if(mod(Sz,2)==0 && Sz>1)
%     M=M + linspace(-pi/2,pi/2,numel(Idxs)).'*ones(1,Sz);
end

Out=exp(1i*M);