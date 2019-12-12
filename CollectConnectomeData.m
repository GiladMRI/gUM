% A=loadniidata('/media/a/DATA1/Downloads/aspera/144832/unprocessed/3T/T1w_MPR1/144832_3T_T1w_MPR1.nii.gz');

% /media/a/DATA1/Downloads/aspera/144832/unprocessed/3T/T1w_MPR2 144832_3T_T1w_MPR2.nii.gz
%%
BaseHCPP='/media/a/DATA/Downloads/aspera/';
DHCP=dir(BaseHCPP);
DHCP=DHCP([DHCP.isdir]);
DHCP=DHCP(3:end);
DHCPN={DHCP.name}';
DHCPN=DHCPN(strlength(DHCPN)==6);

ContrastNames={'T1w_MPR1','T1w_MPR2','T2w_SPC1','T2w_SPC2'};

OutP='~/HomeA/HCPDataset/';
%%
for i=1:numel(DHCPN)
    disp(i);
    CurP=[BaseHCPP DHCPN{i} '/unprocessed/3T/'];
    DP=dir(CurP);
    DP=DP([DP.isdir]);
    DP=DP(3:end);
    DPN={DP.name}';
    for j=1:numel(DPN)
        CurFN=[BaseHCPP DHCPN{i} '/unprocessed/3T/' DPN{j} filesep DHCPN{i} '_3T_' DPN{j} '.nii.gz'];
        A=loadniidata(CurFN);
        for a=1:3
            AP=permute(A,[a setdiff(1:3,a)]);
            for k=1:3
                CurIdx=randi(size(AP,1));
                CurI=squeeze(AP(CurIdx,:,:));
                OutFN=[num2str(i,'%04d') '_' num2str(j,'%01d') '_' num2str(a,'%01d') '_' num2str(k,'%01d') '.mat'];
                CurIc=int32(crop(CurI,[256 256]));
                save([OutP OutFN],'CurIc');
            end
        end
    end
end
%%
