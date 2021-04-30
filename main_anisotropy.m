%% Section 0: Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('/Users/thomas/Documents/MATLAB/added_codes/tiffread2'))
addpath(genpath('/Users/thomas/Documents/MATLAB/Universal_Code_TvZ'))
pathname=uigetdir; pathname=[pathname '/']; cd(pathname)
%% Section 1: Variable and image loading and conversion %%%%%%%%%%%%%%%%%%%
ws=10;%width of the measurement ROIs
bin=3;%data binning pixel-anisotropy determination (should be an odd number)
blur=3;%anisotropy image Gaussian blur for visualisation purposes

temp_PA=tiffread2b([pathname 'PA.tif'],1,40); 
for i=1:length(temp_PA), PE(i,:,:)=temp_PA(i).data; end
temp_PE=tiffread2b([pathname 'PE.tif'],1,40);
for i=1:length(temp_PE), PE(i,:,:)=temp_PE(i).data; end

cam.mode=0; cam.serialNo='0000'; cam.tInt=1000;%camera settings for more info
cam.ROI=[0 0 size(temp_PA.data,1) size(temp_PA.data,2)];%see PhotoConvertIm.m
[PE, gain_PA, offset_PA] = PhotoConvertIm (PE, cam);
[PE, gain_PE, offset_PE] = PhotoConvertIm (PE, cam);

PA_BG=tiffread2([pathname 'BG/BG_PA.tif']);PA_BG=double(PA_BG.data);
PE_BG=tiffread2([pathname 'BG/BG_PE.tif']);PE_BG=double(PE_BG.data);
%BE AWARE TO REMOVE gain and offset IF CAMERA IS NOT EMCCD
BG(1,:,:) = PhotoConvertIm (PA_BG, cam, gain_PA, offset_PA);
BG(2,:,:) = PhotoConvertIm (PE_BG, cam, gain_PE, offset_PE);
PA_GF=tiffread2([pathname 'BG/GF_PA.tif']);PA_GF=double(PA_GF.data);
PE_GF=tiffread2([pathname 'BG/GF_PE.tif']);PE_GF=double(PE_GF.data);
%BE AWARE TO REMOVE gain and offset IF CAMERA IS NOT EMCCD
GF(1,:,:) = PhotoConvertIm (PA_GF, cam, gain_PA, offset_PA);
GF(2,:,:) = PhotoConvertIm (PE_GF, cam, gain_PE, offset_PE);

for i=1:size(PE,1)
    tim1(:,:)=PA(i,:,:);tim2(:,:)=PE(i,:,:);
	PA(i,:,:)=ImageCorrection(tim1,PA_BG);
    [timcorr tFORM(i)] = DualCh_align (tim1, tim2, 'affine', 1);
    [PE(i,:,:) GFactor]=ImageCorrection(tim2,BG,'An',GF, tFORM(i));
    clear tim1 tim2 timcorr
end

PA(find(PA<0))=0;PA(find(PA==Inf))=0;PA(isnan(PA))=0;
PE(find(PE<0))=0;PE(find(PE==Inf))=0;PE(isnan(PE))=0;
clear temp_PA temp_PE PA_BG PE_BG BG PA_GF PE_GF GF
%% Section 2: Indicating the square ROIs on cells to be analysed per image 
for j=1:size(PE,1)
    
   	chPA(:,:)=PA(j,:,:);
    chPE(:,:)=PE(j,:,:);
 	
    No=size(ROI,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(9)
    set(gcf, 'Position',  [400 50 950 900])
    imagesc(chPA+2*chPE);
        title(['Image ' num2str(j) ' of ' num2str(size(PE,1)) ' images: PLEASE DRAW BACKGROUND'])
 
roi = imfreehand%%%Draw a single region that will become the image-associated background

BW=createMask(roi);
BG(1)=sum(sum(BW.*chPA))/sum(sum(BW));
BG(2)=sum(sum(BW.*chPE))/sum(sum(BW));
delete(roi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fin=0;
    while fin==0
       [ROI, fin] = ROIselect(chPA+chPE, ws, 'square', ROI);
    end
        
    for i=No+1:size(ROI,2)
        ROI(i).BG=BG;
        ROI(i).frame=j;
        ROI(i).PA=chPA(ROI(i).area(3):ROI(i).area(4),ROI(i).area(1):ROI(i).area(2));
        ROI(i).PE=chPE(ROI(i).area(3):ROI(i).area(4),ROI(i).area(1):ROI(i).area(2));
    end
    
    im(j).tot=chPA+2*chPE;
    im(j).an=(chPA-chPE)./im(j).tot;    
%% binning the images%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fun = @(x) sum(x(:));
    chPA = nlfilter(chPA,[bin bin],fun); 
        chPA=chPA(1+((bin-1)/2):bin:end-((bin-1)/2),1+((bin-1)/2):bin:end-((bin-1)/2));
    chPE = nlfilter(chPE,[bin bin],fun); 
        chPE=chPE(1+((bin-1)/2):bin:end-((bin-1)/2),1+((bin-1)/2):bin:end-((bin-1)/2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    im(j).BinTot=chPA+2*chPE;
    im(j).BinAn=(chPA-chPE)./im(j).BinTot;
    
    clear No BG BW chPA chPE roi fin
    close Figure 9
    
end
cd(pathname)
save('ROIs.mat')
%% Section 3: Analyzing all the ROIs and extracting relevant data %%%%%%%%
j=1;
for i=1:length(ROI)  
    if length(ROI(i).PA)>0    
    ROI(i).totInt=sum(sum(ROI(i).PA))+2*sum(sum(ROI(i).PE)); 
    ROI(i).stdInt=std2(ROI(i).PA+2*ROI(i).PE);
    ROI(i).meanInt=mean2(ROI(i).PA+2*ROI(i).PE); 
    ROI(i).CountRate=(1000*ROI(i).meanInt)/cam.tInt;
    ROI(i).AN=(sum(sum(ROI(i).PA))-sum(sum(ROI(i).PE)))/ROI(i).totInt;
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%DATA ON BACKGROUND%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=1:length(bleb)
%   BG(i,:)=Bleb(i).BG;
%end
%BGratio=unique(BG(:,1)./BG(:,2));
%BGtotal=unique(BG(:,1)+BG(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
cd(pathname)
save('Analysed_final.mat')
%% Section 4: Data representation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





