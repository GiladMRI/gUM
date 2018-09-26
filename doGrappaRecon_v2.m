function [unfolded_data weights ] = doGrappaRecon_v2(varargin )
% function [unfolded_data weights ] = doGrappaRecon_v2(folded_data, ref_data, acc_fac, kernelsize )      if the kernel has not been fitted yet
% or 
% function [unfolded_data weights ] = doGrappaRecon_v2(folded_data, weights, acc_fac,  kernelsize )      if the kernel is already known
%
% this function performs a box standard GRAPPA reconstruction;  somewhat adapted form Mark Griswohld's OpenGrappa.m 
% INPUT: 
% folded_data    - k-space in format [COL LIN CHA .... SLC .... ]  SLC *must* be in DIM=10
% kernelsize     - size of grappa kernel [ x y ] (x odd and y even numbers)
% acc_fac         - acceleration factor
% ref_data       - ACS reference data as [COL LIN CHA SLC]
% weights        - grappa weights (OPTIONAL -- if they are not passed, they will be fitted) 
% OUTPUT:
% unfolded_data  - reconstructed data (LIN dimension will be acc_fac*LIN of the folded data)
% weights{slice} - Grappa weights for each slice for use in a later reconstruction 
%
% NB for the estimation of the GRAPPA weights, only the *first* value of all additional non-singletom dimensions 
% in the reference data will be considered, e.g. first echo, first set, first phase etc etc...
%
% Benedikt Poser UH March 2013 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we expect the Input to be arranged like this
% DIM_COL = 1;
% DIM_LIN = 2;
% DIM_CHA = 3;
% DIM_SET = 4;
% DIM_ECO = 5;
% DIM_PHS = 6;
% DIM_REP = 7;
% DIM_SEG = 8;
% DIM_PAR = 9;
% DIM_SLC = 10;
% DIM_IDA = 11;
% DIM_IDB = 12;
% DIM_IDC = 13;
% DIM_IDD = 14;
% DIM_IDE = 15;
% DIM_AVG = 16;
% MAX_DIM = 16; 

folded_data = varargin{1};
if prod(size(varargin{2})) > size (folded_data, 10)           % kernel NOT supplied --> doGrappaRecon_v2(folded_data, ref_data, acc_fac, kernelsize )
    ref_data    = varargin{2};
    acc_fac     = varargin{3};
    kernelsize  = varargin{4};
    fit_kernel  = true;
else                                                          % kernel IS supplied --> doGrappaRecon_v2(folded_data, acc_fac, kernelsize )    
    weights     = varargin{2};
    acc_fac     = varargin{3};
    kernelsize  = varargin{4};
    fit_kernel  = false;
end

nCol = size(folded_data,1);
nLin = size(folded_data,2)*acc_fac ;
nCha = size(folded_data,3);
nSlc = size(folded_data,10);

grx = kernelsize(1); 
gry = kernelsize(2);
dx = floor(grx/2);
dy = (gry/2-1)*acc_fac;

%% if NO weights are passed, we calculate them first 
if fit_kernel

    % in case the refdata has more dimensions used that we want here, we only consider the *first* ECO, PHS, SET AVG etc etc ....
    ref_data = squeeze(ref_data(:,:,:,1,1,1,1,1,1,:,1,1,1,1,1,1));
   
    % in this code we need [ CHA LIN COL SLC] order
    ref_data = permute(ref_data, [ 3 2 1 4]); 
    
    nLin_ref = size(ref_data,2);
    nCol_ref = size(ref_data,3);
       
    % loop over slices to determine the weights    
    for slc=1:nSlc
        ref_data_slc = ref_data(:,:,:,slc);
        % very embarassing slow loop to gather all the source and target points for the fit
        source = zeros(nCha*gry*grx,(nLin_ref-(gry-1)*acc_fac)*(nCol_ref-(grx-1)));
        target = zeros(nCha*(acc_fac-1),(nLin_ref-(gry-1)*acc_fac)*(nCol_ref-(grx-1)));
        ii=0;                                                           
        for x=dx+1:nCol_ref-dx,
            for y=1:nLin_ref-(gry-1)*acc_fac,
                ii=ii+1;        
                source(:,ii) = reshape(ref_data_slc(:,y:acc_fac:y+(gry-1)*acc_fac,x-dx:x+dx),nCha*gry*grx,1);
                target(:,ii) = reshape(ref_data_slc(:,y+1+dy:y+dy+acc_fac-1,x),nCha*(acc_fac-1),1);
            end
        end
        % finally fit source points to target points
        weights{slc}=target/source;   
    end
    
end

%% else the kernel has been passed so we dive straight into the recon part 
% permute to [ CHA LIN COL ...... SLC ..... ] order
folded_data = permute(folded_data, [ 3 2 1 4 5 6 7 8 9 10 11 12 13 14 15 16]);  %

% create array for reconstructed data:
target_size = size(folded_data); target_size(2) = nLin; %nLin = size(folded_data,2)*acc_fac ;
target_size(numel(target_size)+1:16)=1;
unfolded_data = complex(zeros(target_size));
for slc=1:nSlc
    for SetCount = 1:target_size(4)
        for EcoCount = 1:target_size(5)
            for PhsCount = 1:target_size(6)
                for RepCount = 1:target_size(7)
                    for SegCount = 1:target_size(8)
                        for ParCount = 1:target_size(9) % not that PAR makes any sense...
                            for IdaCount = 1:target_size(11)
                                for IdbCount = 1:target_size(12)
                                    for IdcCount = 1:target_size(13)
                                        for IddCount = 1:target_size(14)
                                            for IdeCount = 1:target_size(15)
                                                for AvgCount = 1:target_size(16)
                                                    
                tmp_unfolded = complex(zeros(nCha , nLin+2*dy+1 , nCol+2*dx));                          
                tmp_unfolded(: ,(dy+1:acc_fac:nLin+dy) , dx+1:end-dx) ... 
                    = folded_data(:,:,:,SetCount,EcoCount,PhsCount,RepCount,SegCount,ParCount,slc,IdaCount,IdbCount,IdcCount,IddCount,IdeCount,AvgCount)  ;
                for x = dx+1:nCol+dx
                    for y = (1:acc_fac:nLin)
                        % take source points and apply weights to them to fill targets
                        source = reshape(tmp_unfolded(: , y:acc_fac:(y+(gry-1)*acc_fac) , x-dx:x+dx) , [nCha*gry*grx , 1]);
                        tmp_unfolded(:,y+dy+1:y+dy+acc_fac-1,x) = reshape(weights{slc}*source,[nCha (acc_fac-1)]); 
                    end
                end
                % finally: cut out the "good" data that we want to use, and
                % asign to output array with all slices and other dimensions
                 unfolded_data(:,:,:,SetCount,EcoCount,PhsCount,RepCount,SegCount,ParCount,slc,IdaCount,IdbCount,IdcCount,IddCount,IdeCount,AvgCount) ...
                     = tmp_unfolded(:,dy+1:nLin+dy,dx+1:nCol+dx);
                                                                         
                                                    
                                                end
                                            end 
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end 
end %slc

% very finally: permute back to the input/output format [COL LIN CHA ..... SLC....]
unfolded_data = permute (unfolded_data, [3 2 1 4 5 6 7 8 9 10 11 12 13 14 15 16]); 

 

  