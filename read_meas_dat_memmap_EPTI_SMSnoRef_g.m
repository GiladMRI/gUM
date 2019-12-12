
%% 29 October 2009 
% Kawin Setsompop

%% 6 Nov 2009
% shrink the patrefscan_phascor matrix to save room when saving data (line ~ 55)
% fix hard coding of number of phase correct lines assumption (3)
% Kawin Setsompop

% Based on code written by Jennifer McNab and Thomas Witzel on memmap method of
% reading in data quickly 
% read first repetition using read_meas_dat from Jon Polimeni and use memmap
% method to read the rest. 

% six important things: prot,evp,data,data_phascor1d,patrefscan,patrefscan_phascor

%IMPORTANT NOTE: if use Jon's recon chain then remove the deinterleaving
%stuff. 

% the data is arranged in a coil reordered way, except for 7T data which
% will be a bit weird......

%% 6 July 2011
% automatically detect if this file has been read before (i.e. no need for
% the ReadFistRepNeeded_Flag)

%% 19 Jan 2012
% change the way memmap work by using pointer associate with data_begin rather than end of file
% this will allow for reconstuct of incomplete files 
%
% make it work for SMS data 
% make it work for both VB17 and VD11
% NOTE: ReorderChannels option is not included here so will not work with old data that have already got .mat saved of first rep!!!

% using Jon's new read_meas_dat that arrange things in cell array (much faster read)

%% 18 Oct 2016
%
% make the code work for multi-echo data
%% 2018
% make it work for EPTI 
% Fuyixue
% Zijing

function [meas] = read_meas_dat_memmap_EPTI_GESE_SMSnoRef_g(filename,nRepToRead,BeginRep,SMSdata,ascendingSlice_acq,pf_echo,MB_factor)

if nargin == 4
    ascendingSlice_acq = 0; 
end

save_fname = [filename(1:end-4) '_FirstRep.dat'];

if exist([save_fname(1:end-4) '_Raw.mat'],'file') ~= 2
    opt.ReturnStruct=1;
    opt.ReadMultipleRepetitions = 0;
    opt.SqueezeChannels = 1;
    opt.ReturnCellArray = 1; 
    
    % grab meas.data as a cell which make things a lot faster but some cells can be empty due to PAT and SMS

%     opt.SMSRefscan = 1;

    opt.SMSRefscan = 0;
    
    disp('read & save First Rep')
    tic
     
    %readin first rep and save prot,evp,(patrefscan,patrefscan_phascor)
    meas_first = read_meas_dat(filename,opt);
    
    disp('Hack for VD for now!!!*********************************************')
    meas_first.evp.NChaMeas = size(meas_first.data,3);
    disp('*********************************************')
    if isempty(meas_first.evp.RawCol)  % ICE is off when scaning
        meas_first.evp.NLinMeas=size(meas_first.data,2);
        meas_first.evp.NAFLin=meas_first.prot.lAccelFactPE;
%         temp=find(sum(sum(meas_first.data(:,:,1,1,1,1,1,:,1,1),1),8) ~= 0);
%         meas_first.evp.NFirstLin=temp(1);
    %     meas_first.evp.NLinmeas_first 146
    %     meas_first.evp.NFirstRefLin 
        meas_first.evp.NRefLin = meas_first.prot.lRefLinesPE;
        meas_first.evp.NAFPar=1;
        meas_first.evp.NFirstPar=1;
        meas_first.evp.NParMeas=1;
        meas_first.evp.NFirstRefPar=0;
        meas_first.evp.NRefPar=1;
    end
    
    if isfield(meas_first, 'patrefscan') 
        [meas_first.evp,meas_first.prot,meas_first.patrefscan] = SiemensIPATConventionWorkAround(meas_first.evp,meas_first.prot,meas_first.patrefscan);
        % not sure how data is organize esp. when the number of PE lines is not divisible by the acceleration factor
        % Sol: use Jon's convention by looking at non-zero values in his data
        
        % use this if meas.data is a matrix (In Jon's new fast read_meas_dat meas.data is a cell)
        
        if iscell(meas_first.data) == 0 
            meas_first.evp.SegOneLines = find(sum(meas_first.data(:,:,1,1,1,1,1,1,1,1),1) ~= 0);
            meas_first.evp.SegTwoLines = find(sum(meas_first.data(:,:,1,1,1,1,1,2,1,1),1) ~= 0);
        else
        
        % use this if meas.data is a cell
        if SMSdata == 1 % find non empty slice and grab the first SlcGroup slc position
            [SlcMask] = SlcMaskGenerator(meas_first);
            index = find(SlcMask == 1);
            FirstIndex = index(1);           
            %SlcsPerGroup = Nslices/NslicesEX;
            SlcsPerGroup = sum(SlcMask);
            NslicesEX = MB_factor;
            meas_first.prot.SliceSep = sqrt((meas_first.prot.sSliceArray(1).sPosition_dSag- meas_first.prot.sSliceArray(1+SlcsPerGroup).sPosition_dSag)^2+...
            (meas_first.prot.sSliceArray(1).sPosition_dCor - meas_first.prot.sSliceArray(1+SlcsPerGroup).sPosition_dCor)^2+...
            (meas_first.prot.sSliceArray(1).sPosition_dTra - meas_first.prot.sSliceArray(1+SlcsPerGroup).sPosition_dTra)^2);   % distance between two slice

            meas_first.prot.sSliceArray = meas_first.prot.sSliceArray(1:SlcsPerGroup); % CONVENTION: slices position of the top slice group!!!
        else
            FirstIndex = 1;
        end
        b = meas_first.data(1,:,1,1,1,1,1,1,1,FirstIndex);
        b = squeeze(b);
        SegOneLinesMask = length(b);
        for count = 1:length(b)
            SegOneLinesMask(count) = ~isempty(b{count});
        end
        meas_first.evp.SegOneLines = find(SegOneLinesMask == 1);
        
        c = meas_first.data(1,:,1,1,1,1,1,2,1,FirstIndex);
        c = squeeze(c);
        SegTwoLinesMask = length(c);
        for count = 1:length(c)
            SegTwoLinesMask(count) = ~isempty(c{count});
        end
        meas_first.evp.SegTwoLines = find(SegTwoLinesMask == 1); 
        end
    else
        meas_first.evp.SegOneLines = 1:2:meas_first.evp.NLinMeas;
        meas_first.evp.SegTwoLines = 2:2:meas_first.evp.NLinMeas;
        if SMSdata == 1 % find non empty slice and grab the first SlcGroup slc position
            [SlcMask] = SlcMaskGenerator(meas_first);
            index = find(SlcMask == 1);
            FirstIndex = index(1);           
            SlcsPerGroup = sum(SlcMask);
            NslicesEX = MB_factor;
            meas_first.prot.SliceSep = sqrt((meas_first.prot.sSliceArray(1).sPosition_dSag- meas_first.prot.sSliceArray(1+SlcsPerGroup).sPosition_dSag)^2+...
            (meas_first.prot.sSliceArray(1).sPosition_dCor - meas_first.prot.sSliceArray(1+SlcsPerGroup).sPosition_dCor)^2+...
            (meas_first.prot.sSliceArray(1).sPosition_dTra - meas_first.prot.sSliceArray(1+SlcsPerGroup).sPosition_dTra)^2);   % distance between two slice

            meas_first.prot.sSliceArray = meas_first.prot.sSliceArray(1:SlcsPerGroup); % CONVENTION: slices position of the top slice group!!!
        else
            FirstIndex = 1;
        end
        
    end
    
    if size(meas_first.data_phascor1d,5) > 1 % multi-echo data
       disp('find echo containing Nav for each slice')
       % for slc = 1: size(meas_first.data_phascor1d,10)
       eInd = find(meas_first.data_phascor1d(end/2,1,1,1,:,1,1,1,1,1));
       meas_first.data_phascor1d = meas_first.data_phascor1d(:,:,:,:,eInd,:,:,:,:,:);

        %meas_first.data_phascor1d = meas_first.data_phascor1d(:,:,:,:,1,:,:,:,:,:);    
    end
    if (isfield(meas_first, 'patrefscan')) 
        if size(meas_first.patrefscan_phascor,5) > 1 % multi-echo data
        disp('find echo containing Nav for each slice')
        for slc = 1: size(meas_first.patrefscan_phascor,10)
               eInd = find(meas_first.patrefscan_phascor(end/2,1,1,1,:,1,1,1,1,slc));
               meas_first.patrefscan_phascor(:,:,:,:,1,:,:,:,:,slc) = meas_first.patrefscan_phascor(:,:,:,:,eInd,:,:,:,:,slc);
        end
        meas_first.patrefscan_phascor = meas_first.patrefscan_phascor(:,:,:,:,1,:,:,:,:,:);
        
         %disp('assume Nav in last echo')
         %eInd = size(meas_first.patrefscan_phascor,5);
         %meas_first.patrefscan_phascor = meas_first.patrefscan_phascor(:,:,:,:,eInd,:,:,:,:,:);
        
        end
    end
%     if SMSdata == 1
%         if size(meas_first.smsrefscan_phascor,5) > 1 % multi-echo data
%             disp('find echo containing Nav for each slice')
%             for slc = 1: size(meas_first.smsrefscan_phascor,10)
%                eInd = find(meas_first.smsrefscan_phascor(end/2,1,1,1,:,1,1,1,1,slc));
%                meas_first.smsrefscan_phascor(:,:,:,:,1,:,:,:,:,slc) = meas_first.smsrefscan_phascor(:,:,:,:,eInd,:,:,:,:,slc);
%             end
%             meas_first.smsrefscan_phascor = meas_first.smsrefscan_phascor(:,:,:,:,1,:,:,:,:,:);
%            
%         end
%     end
    
    meas_first.evp.PhasCorSegSwap = (sum(meas_first.data_phascor1d(:,1,1,1,1,1,1,1,1,1),1) == 0); % for some data set, Jon's code swap the phase correction segment to match data

    meas_first.evp.NPhaseCorLines = size(meas_first.data_phascor1d,2);
    
    if ascendingSlice_acq == 0
        deinterleave = strcmp(meas_first.prot.ucMultiSliceMode, 'MSM_INTERLEAVED');
    else
        meas_first.prot.ucMultiSliceMode = 'MSM_SEQUENTIAL';
        deinterleave = 2;
    end
    
    if  (isfield(meas_first, 'patrefscan_dpg')) % MKM 
        meas_first.evp.NFirstRefPar = find(sum(meas_first.patrefscan_dpg(end/2,:,1,1,1,1,1,:,1,1),8),1);
        [~, meas_first.patrefscan_dpg] = mrir_array_GRAPPA_prune([], meas_first.patrefscan_dpg , meas_first.evp); % need to modify Jon's code a bit to use this, o.w. use code below
    end
    
    if deinterleave == 1
%         if SMSdata == 1
%             meas_first.smsrefscan = mrir_image_slice_deinterleave(meas_first.smsrefscan);
%             meas_first.smsrefscan_phascor = mrir_image_slice_deinterleave(meas_first.smsrefscan_phascor);
%         end
        if  (isfield(meas_first, 'patrefscan'))
            meas_first.patrefscan = mrir_image_slice_deinterleave(meas_first.patrefscan);
            meas_first.patrefscan_phascor = mrir_image_slice_deinterleave(meas_first.patrefscan_phascor);
        end
         if  (isfield(meas_first, 'patrefscan_dpg'))
            meas_first.patrefscan_dpg = mrir_image_slice_deinterleave(meas_first.patrefscan_dpg);
         end
        if (isfield(meas_first, 'offline'))
            meas_first.offline = mrir_image_slice_deinterleave(meas_first.offline);
        end
    elseif deinterleave == 2 % ascending need to reverse
%         if SMSdata == 1
%             meas_first.smsrefscan = meas_first.smsrefscan(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
%             meas_first.smsrefscan_phascor = meas_first.smsrefscan_phascor(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
%         end
        if  (isfield(meas_first, 'patrefscan'))
            meas_first.patrefscan = meas_first.patrefscan(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
            meas_first.patrefscan_phascor = meas_first.patrefscan_phascor(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
        end
        if  (isfield(meas_first, 'patrefscan_dpg'))
            meas_first.patrefscan_dpg = meas_first.patrefscan_dpg(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
        end
        if (isfield(meas_first, 'offline'))
            meas_first.offline = meas_first.offline(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
        end
    end
    
    meas_first.data = [];
    meas_first.data_phascor1d = []; 
    
    SaveRawData(meas_first,save_fname);
    
    toc
    
else
    disp('Load First Rep Info')
    tic
    meas_first = ReadRawData(save_fname);
    toc
    if ascendingSlice_acq == 0
        deinterleave = strcmp(meas_first.prot.ucMultiSliceMode, 'MSM_INTERLEAVED');
    else
        meas_first.prot.ucMultiSliceMode = 'MSM_SEQUENTIAL';
        deinterleave = 2;
    end  
end

coil_index = 1:meas_first.evp.NChaMeas;

%% extract params from first rep

meas = meas_first;
clear meas_first;
if SMSdata == 1
    if meas.isvd == 1    
        if ~isempty(meas.prot.sWiPMemBlock_adFree)
            %MKM meas.prot.sWiPMemBlock_adFree(1) = meas.prot.sWipMemBlock_adFree(5) % SMS factor
            %MKM meas.prot.sWiPMemBlock_adFree(2) = meas.prot.sWipMemBlock_adFree(3)*meas.prot.lAccelFactPE % FOV shift factor
            meas.prot.sWiPMemBlock_adFree(1) = meas.prot.sWipMemBlock_adFree(3);%*meas.prot.lAccelFactPE % FOV shift factor
            meas.prot.sWiPMemBlock_adFree(2) = meas.prot.sWipMemBlock_adFree(2); % SMS factor
            
        else
            %MKM 
%             disp('cant find SMS and shift paras in .dat')
%             disp('manual entering here.....')
%             keyboard
%             meas.prot.sWiPMemBlock_adFree(1) = 4; % SMS
%             meas.prot.sWiPMemBlock_adFree(2) = 2; % shift in reduced FOV
            [prottemp] = read_meas_prot(filename);
            meas.prot.sWiPMemBlock_adFree(2) = prottemp.sWiPMemBlock_adFree(3); % shift in reduced FOV
            meas.prot.sWiPMemBlock_adFree(1) = prottemp.sWiPMemBlock_adFree(5); % SMS  !!%% original is (2)
            
        end
    else
        %MKM comment out and find manually 
%         meas.prot.sWiPMemBlock_adFree(1) = meas.prot.sWiPMemBlock_adFree(2); % SMS factor
%         meas.prot.sWiPMemBlock_adFree(2) = meas.prot.sWiPMemBlock_adFree(3)*meas.prot.lAccelFactPE; % FOV shift factor
%         meas.prot.sWiPMemBlock_adFree(3) = meas.prot.dThickness/meas.prot.sWiPMemBlock_adFree(1); % slice seperation in mm
        %meas.prot.sWiPMemBlock_adFree
        meas.prot.sWiPMemBlock_adFree(1) = 2; % SMS factor
        meas.prot.sWiPMemBlock_adFree(2) = 2; % FOV shift factor
        meas.prot.sWiPMemBlock_adFree(3) = meas.prot.dThickness/meas.prot.sWiPMemBlock_adFree(1); % slice seperation in mm
    end
    %MKM - hard code in for bay8 scans 
    if meas.prot.sWiPMemBlock_adFree(1) > 10 %
        meas.prot.sWiPMemBlock_adFree(1) = meas.prot.sWiPMemBlock_adFree(2); % SMS factor
        meas.prot.sWiPMemBlock_adFree(2) = meas.prot.sWiPMemBlock_adFree(5); % FOV shift factor
    end
    NslicesEX =  meas.prot.sWiPMemBlock_adFree(1);
else
    NslicesEX = 1;
end

if isempty(meas.evp.NColMeas)==0
    sData = [meas.evp.NColMeas, meas.evp.NLinMeas, meas.evp.NChaMeas,...
              1, meas.evp.NEcoMeas, 1, nRepToRead, 2, 1, meas.evp.NSlcMeas/NslicesEX ];
    sPhaseCor = sData; sPhaseCor(2) = meas.evp.NPhaseCorLines; sPhaseCor(5) = 1;

    nRead = sData(1);
    nPE = sData(2);
    nCoil = sData(3);
    nSlice = sData(10);
    nPhaseCor = sPhaseCor(2);
    nEco = sData(5);

    nPE_Raw = meas.evp.RawLin;
    nRep = meas.evp.RawRep;
    meas.evp.NRepMeas = nRepToRead;
else     
     sData = [meas.prot.lBaseResolution*2, meas.evp.NLinMeas, meas.evp.NChaMeas,...
     1, meas.evp.NEcoMeas, 1, nRepToRead, 2, 1, meas.evp.NSlcMeas/NslicesEX ];
    sPhaseCor = sData; sPhaseCor(2) = meas.evp.NPhaseCorLines;
    sPhaseCor(5)=1; % echo = 1
    if pf_echo==0
        sData(5)=1;
    end
    nRead = sData(1);
    nPE = sData(2);
    nCoil = sData(3);
    nSlice = sData(10);
    nPhaseCor = sPhaseCor(2);
    nEco = sData(5);
    nPE_Raw = size(meas.evp.SegOneLines,2)+size(meas.evp.SegTwoLines,2);
    nRep=nRepToRead;
    meas.evp.NRepMeas = nRepToRead;
end

SegOneLines = meas.evp.SegOneLines;
SegTwoLines = meas.evp.SegTwoLines;
PhasCorSegSwap  = meas.evp.PhasCorSegSwap;


%% calculate parameters needed for mmap and readout and preallocate matrices

meas.data = single(zeros(sData));
meas.data_phascor1d = single(zeros(sPhaseCor));

if pf_echo>0
    Npe_1echo= floor(pf_echo*nPE_Raw);
    Npe_2echo= nPE_Raw;
    linesperrep = round(Npe_2echo+Npe_1echo+nPhaseCor)*nCoil*nSlice;
else
    linesperrep=round(nPE_Raw+nPhaseCor)*nCoil*nSlice ;
end
%% mmap
disp('Mmap')
tic

%keyboard;
%nRep = 5;
nRep=max(nRep,BeginRep);
m = mmap_mdh_noheader_rep_v2(filename,linesperrep,nRep,nRead,nCoil,meas.data_begin,meas.isvd);
% note: meas.isvd is added into read_meas_dat by Kawin to have flag to see if vd or vb
%m = mmap_mdh_noheader_rep(filename,linesperrep,nRep);
disp(['Time: ' num2str(toc) ' s'])


%% Readout

ICE_RAWDATA_SCALE       = 131072.0;  % 64 ^ 3 / 2
K_ICE_AMPL_SCALE_FACTOR = 80 * 20 * ICE_RAWDATA_SCALE / 65536;
if meas.isvd == 1
    mdh_length_float32 = 48; % per PE line 
    mdh_ch_length_float32 = 8; % per chanel in each PE line
else
    mdh_length_float32 = 0;
    mdh_ch_length_float32 = 32;
end


for RepCount = 1:nRepToRead
    disp([ 'Reading in Rep:' num2str(BeginRep+RepCount-1) ])
    tic
    k = m.Data(BeginRep+(RepCount-1)).kdata;
%     k = reshape(k,[],linesperrep)*K_ICE_AMPL_SCALE_FACTOR;
%     k = k(33:2:end,:)+i*k(34:2:end,:); %32 headers
%     k = reshape(k, [nRead, nCoil, nPE_Raw+nPhaseCor, nSlice]);
    k = reshape(k,[],linesperrep/nCoil)*K_ICE_AMPL_SCALE_FACTOR;   
    k = k(1+mdh_length_float32:end,:);
    k = reshape(k,[],linesperrep);    
    k = k(1+mdh_ch_length_float32:2:end,:) + 1i*k(2+mdh_ch_length_float32:2:end,:); %32 headers
    if pf_echo>0
        k = reshape(k, [nRead, nCoil, (Npe_1echo+Npe_2echo)+ nPhaseCor, nSlice]);
    else
        k = reshape(k, [nRead, nCoil,  nPE_Raw + nPhaseCor, nSlice]);
    end

    if PhasCorSegSwap == 0
        meas.data_phascor1d(:,1:2:end,:,:,:,:,RepCount,1,:,:) = permute(k(:,coil_index,1:2:nPhaseCor,:),[1 3 2 4]);
        meas.data_phascor1d(:,2:2:end,:,:,:,:,RepCount,2,:,:) = permute(k(end:-1:1,coil_index,2:2:nPhaseCor,:),[1 3 2 4]); % need to reverse k-line as assume EPI data
    else
        meas.data_phascor1d(:,2:2:end,:,:,:,:,RepCount,1,:,:) = permute(k(:,coil_index,2:2:nPhaseCor,:),[1 3 2 4]);
        meas.data_phascor1d(:,1:2:end,:,:,:,:,RepCount,2,:,:) = permute(k(end:-1:1,coil_index,1:2:nPhaseCor,:),[1 3 2 4]); % need to reverse k-line as assume EPI data
    end
    
    if pf_echo>0
        k0=zeros(nRead,nCoil,Npe_2echo,nEco,nSlice);
        k0(:,:,1:Npe_1echo,1,:) = reshape(k(:,:,nPhaseCor+1:nPhaseCor+Npe_1echo,:),nRead,nCoil,Npe_1echo,1,nSlice); 
        k0(:,:,:,2,:) = reshape(k(:,:,nPhaseCor+Npe_1echo+1:end,:),nRead,nCoil,Npe_2echo,1,nSlice); 
    else
        k0=zeros(nRead,nCoil,nPE_Raw,nEco,nSlice);
        k0(:,:,1:nPE_Raw,:,:) = reshape(k(:,:,nPhaseCor+1:end,:),nRead,nCoil,nPE_Raw,1,nSlice); 
    end
    k=k0;
    
%    if (isfield(meas, 'offline') && mod(meas.prot.lAccelFactPE,2))
%        meas.data(:,SegOneLines,:,:,:,:,RepCount,1,:,:) = permute(k(end:-1:1,coil_index,1:2:end,:,:),[1 3 2 4 5]);
%        meas.data(:,SegTwoLines,:,:,:,:,RepCount,2,:,:) = permute(k(:,coil_index,2:2:end,:,:),[1 3 2 4 5]); % need to reverse k-line as assume EPI data
%    else
        meas.data(:,SegOneLines,:,:,:,:,RepCount,1,:,:) = permute(k(:,coil_index,1:2:end,:,:),[1 3 2 4 5]);
        meas.data(:,SegTwoLines,:,:,:,:,RepCount,2,:,:) = permute(k(end:-1:1,coil_index,2:2:end,:,:),[1 3 2 4 5]); % need to reverse k-line as assume EPI data
%   end
    toc
end

disp(' ')
disp(' ')


%if ( deinterleave ),
%    meas.data = mrir_image_slice_deinterleave(meas.data);
%    meas.data_phascor1d = mrir_image_slice_deinterleave(meas.data_phascor1d);
%end
%deinterleave =0; % MKM - test 

if deinterleave == 1
    meas.data = mrir_image_slice_deinterleave(meas.data);
    meas.data_phascor1d = mrir_image_slice_deinterleave(meas.data_phascor1d);
elseif deinterleave == 2 % ascending need to reverse
    meas.data = meas.data(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
    meas.data_phascor1d = meas.data_phascor1d(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
end

%MKM 
if  (isfield(meas, 'patrefscan'))
   % meas.evp.NFirstRefPar = find(sum(meas.patrefscan(end/2,:,1,1,1,1,1,:,1,1),8),1); %MKM 
    [meas.data] = mrir_array_GRAPPA_prune(meas.data, [] , meas.evp);
end
clear m
    
function [evp,prot,patrefscan] = SiemensIPATConventionWorkAround(evp,prot,patrefscan) 
%% Modify MEAS file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% workaround for unusual Siemens convention #1:
if ( evp.NFirstRefLin == 0 ),
    evp.NFirstRefLin = mrir_ice_dimensions(patrefscan, 'lin') - evp.NRefLin + 1;
end;

% workaround for unusual Siemens convention #2:

% (possibly Siemens fills in with GRAPPA fewer lines than are in FFT, so
% last lines are effectively zero-padded; this could throw off SNR
% calculations, so by overriding this we force "mrir_epi_GRAPPA" to fill
% in same number of k-space lines as there are image lines.)
if ( ~isempty(evp.NAFLin) && (evp.NAFLin == 1) && (prot.ucPhasePartialFourier == 1) && (evp.NLinMeas < evp.NImageLins) ),
    jnotify
    keyboard
    evp.NLinMeas = evp.NImageLins;
end;

function [SlcMask] = SlcMaskGenerator(meas_first)

% figure out SMS factor and slc location

if  isfield(meas_first, 'patrefscan')
    R_inplane = meas_first.prot.lAccelFactPE;
else
    R_inplane = 1;
end

PElineSearchCount = 1; 
while PElineSearchCount <= R_inplane+1
    a = meas_first.data(1,PElineSearchCount,1,1,1,1,1,1,1,:);
    a = squeeze(a);
    SlcMask = length(a);
    for count = 1:length(a)
        SlcMask(count) = ~isempty(a{count});
    end
    if sum(SlcMask) ~= 0
        break; % found data!!
    else
        if (PElineSearchCount == R_inplane+1)
            disp('error!! cant find the data')
            keyboard
        end
        PElineSearchCount = PElineSearchCount+ 1; % keep looking!!
    end
end

%SlcMask = find(sum(meas_first.data(end/2,:,1,1,1,1,1,2,1,:),2) ~= 0)

%SMS_factor = size(meas_first.smsrefscan,10)/sum(SlcMask);
    
    
    
    
    
     
  
  
  
  

