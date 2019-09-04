function gsafe_plot2(g_mTm,GradDwellTime_ms,Gsystems,GradPerm)
if(nargin<3)
    Gsystems=[];
end
if(nargin<4)
    GradPerm=[];
end
if(isempty(Gsystems))
    Gsystems={safe_params_prisma() safe_params_ma7T() safe_params_trio() safe_params_skyra()};
end
if(isstruct(Gsystems))
    Gsystems={Gsystems};
end
if(isempty(GradPerm))
    GradPerm=1:3;
end

% figure;

gwf=[real(g_mTm) imag(g_mTm) imag(g_mTm)*0]/1000;

gwf=gwf(:,GradPerm); % For Sagittal, Coronal

dt=GradDwellTime_ms/1000;
% Predict PNS levels
pnss=cell(1,numel(Gsystems));
Ttls=cell(1,numel(Gsystems));
for i=1:numel(Gsystems)
    pnss{i}=safe_gwf_to_pns(gwf, dt, Gsystems{i});
    Ttls{i}=Gsystems{i}.name;
end
% Plot some results
% gsafe_plot([gmat2cell(pnsPrisma(:,1:2),2) gmat2cell(pnsma7T(:,1:2),2)], dt);

% if nargin < 2
%     ttot    = 1; % au
%     xlabstr = 'Time [a.u.]';
% else
%     ttot = size(pns,1) * dt * 1000; % ms
    xlabstr = 'Time [ms]';
% end

% h = plot(t,pns); 
CLRs='rgbcmrgbcmkrgbcm';
Mrkrs={'-','--',':'};
k=1;
clear h Lgnds
for j=1:numel(pnss)
    ttot = size(pnss{j},1) * dt * 1000; % ms
    t  = linspace(0, ttot, size(pnss{j},1));
    for i=1:3
        h(k)=plot(t,pnss{j}(:,i),[CLRs(i) Mrkrs{j}]); hold on;
        Lgnds{k}=[Ttls{j} ' ' num2str(max(pnss{j}(:,i)),'%0.0f') '%'];
        k=k+1;
    end
    JointPNS=sqrt(sum(pnss{j}.^2,2));
    h(k)=plot(t,JointPNS,['k' Mrkrs{j}]); hold on;
    Lgnds{k}=[Ttls{j} ' ' num2str(max(JointPNS),'%0.0f') '%'];
    k=k+1;
end

hold on; grid on;

ylim([0 150])
xlim([min(t) max(t)])

title('Predicted PNS');

xlabel (xlabstr)
ylabel ('Relative stimulation [%]')

% h(2) = plot([0 max(t)], [1 1] * max(pns(:)), 'k--');
legend(h,Lgnds);
% C={};
% for i=1:numel(pns)
%     C{i}=[num2str(i) ' (' num2str(max(pns{i}), '%0.0f') '%)'];
% end
% h(3) = legend(C{:},'location', 'best');
% h(3) = legend(...
%     ['X (' num2str(max(pns(:,1)), '%0.0f') '%)'], ...
%     ['Y (' num2str(max(pns(:,2)), '%0.0f') '%)'], ...
%     ['Z (' num2str(max(pns(:,3)), '%0.0f') '%)'], ...
%     'location', 'best');