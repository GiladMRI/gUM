function Out=GenerateRandomSinPhase(N,LFac,QFac,SFac,RandV)
if(numel(N)==1)
    N=[N N];
end
if(nargin<5)
    RandV=rand(11,1);
end
if(nargin<3)
    QFac=1;
end
if(nargin<2)
    LFac=5;
end
if(nargin<4)
    SFac=1;
end
Linx=linspace(-pi,pi,N(1));
Liny=linspace(-pi,pi,N(2));

[X,Y]=ndgrid(Linx,Liny);

AL=(RandV(1)-0.5)*LFac*X+(RandV(2)-0.5)*LFac*Y+(RandV(3)-0.5)*QFac*(X.^2)+(RandV(4)-0.5)*QFac*(Y.^2);
BL=(RandV(5)-0.5)*LFac*X+(RandV(6)-0.5)*LFac*Y+(RandV(7)-0.5)*QFac*(X.^2)+(RandV(8)-0.5)*QFac*(Y.^2);

PX=(RandV(9)-0.5)*sin(AL)+(RandV(10)-0.5)*sin(BL);

DCPhase=RandV(11)*2*pi-pi;
PX=PX*2*pi*SFac+DCPhase;

Out=exp(1i*PX);

%%
% # LFac = 5
% 	# QFac = 0.1
% 	# nx, ny = (3, 2)
% 	Linx = tf.linspace(-np.pi, np.pi, nx)
% 	Liny = tf.linspace(-np.pi, np.pi, ny)
% 	X, Y = tf.meshgrid(Linx, Liny,indexing='ij')
% 
% 	Rnd=tf.random_uniform([11])
% 
% 	AL=(Rnd[0]-0.5)*LFac*X+(Rnd[1]-0.5)*LFac*Y+(Rnd[2]-0.5)*QFac*( tf.pow(X,2)  )+(Rnd[3]-0.5)*QFac*(  tf.pow(Y,2)  );
% 	BL=(Rnd[4]-0.5)*LFac*X+(Rnd[5]-0.5)*LFac*Y+(Rnd[6]-0.5)*QFac*( tf.pow(X,2)  )+(Rnd[7]-0.5)*QFac*(  tf.pow(Y,2) );
% 	PX=(Rnd[8]-0.5)*tf.sin(AL)+(Rnd[9]-0.5)*tf.sin(BL);
% 	DCPhase=Rnd[10]*2*np.pi-np.pi;
% 	PX=PX*2*SFac*np.pi+DCPhase;
% 	Out=tf.exp(tf.complex(PX*0, PX));
% 	return Out