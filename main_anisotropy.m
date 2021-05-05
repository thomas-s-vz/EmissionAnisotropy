%main_anisotropy has been written by Thomas S van Zanten last update 5th May 2021
%% Section 0: Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('path where all the general functions are located'))
addpath(genpath('path where all the specific functions are located'))
pathname=uigetdir; pathname=[pathname '/']; cd(pathname)
%% Section 1: Variable and image loading and conversion %%%%%%%%%%%%%%%%%%%
ws=10;%width of the measurement ROIs
bin=3;%data binning pixel-anisotropy determination (should be an odd number)
blur=3;%anisotropy image Gaussian blur for visualisation purposes

temp_PA=tiffread2b([pathname 'PA.tif'],1,40); 
for i=1:length(temp_PA), PA(i,:,:)=temp_PA(i).data; end
temp_PE=tiffread2b([pathname 'PE.tif'],1,40);
for i=1:length(temp_PE), PE(i,:,:)=temp_PE(i).data; end

cam.mode=0; cam.serialNo='0000'; cam.tInt=2000;%camera settings for more info
cam.ROI=[0 0 size(PA,2) size(PA,3)];%see PhotoConvertIm.m
[PA, gain_PA, offset_PA] = PhotoConvertIm (PA, cam);
[PE, gain_PE, offset_PE] = PhotoConvertIm (PE, cam);

PA_BG=tiffread2([pathname 'BG/BG_PA.tif']);PA_BG=double(PA_BG.data);
PE_BG=tiffread2([pathname 'BG/BG_PE.tif']);PE_BG=double(PE_BG.data);
%BE AWARE TO REMOVE gain and offset IF CAMERA IS NOT EMCCD
PA_BG = PhotoConvertIm (PA_BG, cam, gain_PA, offset_PA);
BG(1,:,:) = PA_BG;
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
    [PE(i,:,:) GFactor]=ImageCorrection(tim2, BG, 'Ani', GF, tFORM(i));
    clear tim1 tim2 timcorr
end

PA(find(PA<0))=0;PA(find(PA==Inf))=0;PA(isnan(PA))=0;
PE(find(PE<0))=0;PE(find(PE==Inf))=0;PE(isnan(PE))=0;
clear temp_PA temp_PE PA_BG PE_BG BG PA_GF PE_GF GF
cd(pathname)
save('corrected_images.mat')
%% Section 2: Indicating the square ROIs on cells to be analysed per image 
c=1; ROI=[];
for j=1:size(PE,1)
    
   	chPA(:,:)=PA(j,:,:);
    chPE(:,:)=PE(j,:,:);
 	
    No=size(ROI,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(9)
    set(gcf, 'Position',  [400 50 950 900])
    imagesc(chPA+2*chPE);
        title(['Image ' num2str(j) ' of ' num2str(size(PE,1)) ' images: PLEASE DRAW BACKGROUND'])
 
imroi = imfreehand%%%Draw a single region that will become the image-associated background
setColor(imroi,'black');
BW=createMask(imroi);
BG(1)=sum(sum(BW.*chPA))/sum(sum(BW));
BG(2)=sum(sum(BW.*chPE))/sum(sum(BW));
cellBW=zeros(size(chPA));
title(['Image ' num2str(j) ' of ' num2str(size(PE,1))...
    ' images: PLEASE DRAW THE CELLS (double-click in empty space when done)'])
    while sum(sum(BW))>100
        imroi = imfreehand%%%Draw the cells in the image
        setColor(imroi,'green');
        BW=createMask(imroi);
        cellBW=cellBW + c*BW; c=c+1;
    end
c=c-1; delete(imroi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fin=0;
    while fin==0
       [ROI, fin] = ROIselect(chPA+chPE, ws, 'square', ROI);
    end
        
    for i=No+1:size(ROI,2)
        ROI(i).BG=BG;
        ROI(i).frame=j;
        ROI(i).cell=cellBW(ROI(i).area(3)+ws/2,ROI(i).area(1)+ws/2);
        ROI(i).PA=chPA(ROI(i).area(3):ROI(i).area(4),ROI(i).area(1):ROI(i).area(2));
        ROI(i).PE=chPE(ROI(i).area(3):ROI(i).area(4),ROI(i).area(1):ROI(i).area(2));
    end
    
    im(j).cells=cellBW;
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
    
    clear No BG BW chPA chPE imroi fin cellBW
    close Figure 9
    
end
clear c j i
cd(pathname)
save('images_ROIs.mat', 'ROI', 'im')
%% Section 3: Analyzing all the ROIs and extracting relevant data %%%%%%%%
j=1;
for i=1:length(ROI)  
    if length(ROI(i).PA)>0    
    ROI(i).totInt=sum(sum(ROI(i).PA))+2*sum(sum(ROI(i).PE)); 
    ROI(i).stdInt=std2(ROI(i).PA+2*ROI(i).PE);
    ROI(i).meanInt=mean2(ROI(i).PA+2*ROI(i).PE); 
    ROI(i).CountRate=(1000*ROI(i).meanInt)/cam.tInt;
    ROI(i).AN=(sum(sum(ROI(i).PA))-sum(sum(ROI(i).PE)))/ROI(i).totInt;
%error in anisotropy from photon statistics error propagation: Lidke et al 2005
    ROI(i).Err=((1-ROI(i).AN)*(1+2*ROI(i).AN)*(1-ROI(i).AN+mean(mean(GFactor))*...
        (1+2*ROI(i).AN)))/(3*ROI(i).totInt);
    IntAn(i,:)=[ROI(i).meanInt ROI(i).AN ROI(i).Err/ROI(i).AN];
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%DATA ON BACKGROUND%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=1:length(ROI)
%   BG(i,:)=ROI(i).BG;
%end
%BGtotal=unique(BG(:,1)+BG(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
cd(pathname)
save('Analysed_final.mat')
%% Section 4: Data representation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=8; figure(i)
set(gcf, 'Position', [0 50 1500 1000]);
% subplot(2,2,1)
% imagesc(im(i).tot), colormap(gca,'gray'), axis off, axis square, caxis([0 0.8*max(max(im(i).tot))]), colorbar
% title('Total Intensity','FontSize',16)
%
% subplot(2,2,2)
% imagesc(im(i).an), colormap(gca,'jet'), axis off, axis square, caxis([0.1 0.3]), colorbar
% title(Anisotropy,'FontSize',16)

subplot(2,2,1)
imagesc(im(i).BinTot), colormap(gca,'gray'), axis off, axis square, caxis([0 0.8*max(max(im(i).BinTot))]), colorbar
title(['Total Intensity from bin ' num2str(bin)],'FontSize',16)

subplot(2,2,2)
imagesc(im(i).BinAn), colormap(gca,'jet'), axis off, axis square, caxis([0.1 0.3]), colorbar
title(['Anisotropy from bin ' num2str(bin)],'FontSize',16)

subplot(2,2,3:4)
plot(IntAn(:,1), IntAn(:,2), 'bo')
axis([0.8*min(IntAn(:,1)) 1.05*max(IntAn(:,1)) 0.1 0.3])
xlabel('Mean intensity (Photons)','FontSize',12), ylabel('Anisotropy','FontSize',12)
title('Total Intensity versus Anisotropy','FontSize',16)

%subplot(2,2,4)
%plot(ws*ws*IntAn(:,1), IntAn(:,3), 'ro')
%axis([0.8*min(IntAn(:,1)) 1.2*max(IntAn(:,1)) 0.8*min(IntAn(:,3)) 1.05*max(IntAn(:,3))])
%title('Associated fractional error in anisotropy calculation','FontSize',16)
