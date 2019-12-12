function Q=CombineRIFlds(Q)

Flds=fieldnames(Q);
SFlds=sort(Flds);

for i=1:numel(SFlds)
    if(strcmp(SFlds{i}(end-2:end),'I_0'))
        Q.([SFlds{i}(1:end-3) 'C'])=double(Q.([SFlds{i}(1:end-3) 'R_0'])+1i*Q.([SFlds{i}(1:end-3) 'I_0']));
        Q=rmfield(Q,[SFlds{i}(1:end-3) 'R_0']);
        Q=rmfield(Q,[SFlds{i}(1:end-3) 'I_0']);
    end
    
    if(SFlds{i}(end)=='I')
        Q.([SFlds{i}(1:end-1) 'C'])=double(Q.([SFlds{i}(1:end-1) 'R'])+1i*Q.([SFlds{i}(1:end-1) 'I']));
        Q=rmfield(Q,[SFlds{i}(1:end-1) 'R']);
        Q=rmfield(Q,[SFlds{i}(1:end-1) 'I']);
    end
end

