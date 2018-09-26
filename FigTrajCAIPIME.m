Traj3=BARTTrajMS;

CAIPIkI = interp1(1:size(kTraj,1),cumsum(CAIPIVec),1:1e5/spBW:(size(kTraj,1)-0.01));

Traj3(3,:)=CAIPIkI;
figure;plot3(Traj3(1,:),Traj3(2,:),Traj3(3,:))
B=CAIPIkI<(max(CAIPIkI)+min(CAIPIkI)/2);
Traj3X=Traj3;
Traj3X(:,B)=NaN;
figure;plot(Traj3X(1,:),Traj3X(2,:),'LineWidth',2);hold on
Traj3X=Traj3;
Traj3X(:,~B)=NaN;
plot(Traj3X(1,:),Traj3X(2,:),'LineWidth',2);
setXaxis([-1.1 1.1]*ceil(MaxK));
setYaxis([-1.1 1.1]*ceil(MaxK));
%%
B1=CAIPIkI>(max(CAIPIkI)-0.1);
B2=CAIPIkI<(min(CAIPIkI)+0.1);
Traj3X=Traj3;
Traj3X(:,~B1)=NaN;
figure;plot3(Traj3X(1,:),Traj3X(2,:),Traj3X(3,:),'LineWidth',4)
hold on
Traj3X=Traj3;
Traj3X(:,~B2)=NaN;
plot3(Traj3X(1,:),Traj3X(2,:),Traj3X(3,:),'LineWidth',4)
Traj3X=Traj3;
Traj3X(:,B1 | B2)=NaN;
plot3(Traj3X(1,:),Traj3X(2,:),Traj3X(3,:),'g--','LineWidth',2)
%%


figure;
plot3(1:1281,kTraj(1:1281,1).',kTraj(1:1281,2).','-','LineWidth',2)
hold on
plot3(1281+(1:1281),kTraj(1281+(1:1281),1).',kTraj(1281+(1:1281),2).','-','LineWidth',2)
xlabel('time');
ylabel('k_x');
zlabel('k_y');
% setXaxis([-1.1 1.1]*ceil(MaxK));
