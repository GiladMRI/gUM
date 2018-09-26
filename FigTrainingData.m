% FigTrainingData
P2DFPBase='/media/a/DATA/P2DF/';
D=dir(P2DFPBase);
D=D(4:end);
%%
for i=1:numel(D)
    Di=dir([P2DFPBase D(i).name filesep '*.*']);
    Di=Di([Di.isdir]);
    Ps{i}=[D(i).name filesep Di(4).name];
end
% Ps={'RegridTry3C2_7TS_P2DF_3.00__2018-07-17_09-23-45_train','RegridTry3C2_7TS_P2DF__2018-07-16_19-19-28_train',...
%     'RegridTry3C2_7TS_RL_S3__2018-07-16_15-19-07_train','RegridTry3C2_7TS_GL_S3__2018-07-10_19-13-55_train'};
%%
for i=1 %:numel(Ps)
    CurP=[P2DFPBase Ps{i} filesep];
    D=dir([CurP 'batch*.png']);
    X=imread([CurP D(1).name]);X=X(:,:,1);
    Y=X(:,[801:1000 1201:1400]);
end
% Z=CombineDims(Y,[3 2]);
Z=PartitionDim(Y,1,4);
Z=CombineDims(Z,[3 2]);