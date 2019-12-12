function AddTopTtls(Ttls,XSz)
MSz=size(get(get(gca,'Children'),'CData'));
% Ttls={'EPTI 1shot','EPTI 1shot','SKEPTIC 1shot','SKEPTIC 3shot'};
for i=1:numel(Ttls)
    h=text(MSz(2)/XSz(2)*i-MSz(2)/XSz(2)/2,7,Ttls{i},'FontSize',14,'Color',[1 0 0.5],'HorizontalAlignment','center'); % set(h,'Rotation',90);
end

% h=text(7,60,'EPTI-3','FontSize',20,'Color',[1 0 0.5],'HorizontalAlignment','center');set(h,'Rotation',90);
% h=text(7,180,'Spiral-3','FontSize',20,'Color',[1 0 0.5],'HorizontalAlignment','center');set(h,'Rotation',90);
