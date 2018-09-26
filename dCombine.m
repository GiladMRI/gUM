for sli=1:12
    RecTSAllX(:,:,sli,:)=RecTSAll(:,:,1,sli,:);
    RecTSAllX(:,:,sli+12,:)=RecTSAll(:,:,2,sli,:);
end
