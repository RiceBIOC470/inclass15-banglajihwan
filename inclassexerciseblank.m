%GB coments
Step1 100
Step2 100
Step3 100 
Step4 100
Step5 100
Step6 100 
Step7 100 
Step8 0 No code or discussion 
Overall 88


%% step 1: write a few lines of code or use FIJI to separately save the
% nuclear channel of the image Colony1.tif for segmentation in Ilastik
file1 = '48hColony1_DAPI.tif' 
reader = bfGetReader(file1); 
reader.getSizeC
chan = 1;
zplane = 1;
time = 1;
iplane1 = reader.getIndex(chan-1, zplane-1, time-1) +1 
img1 = bfGetPlane(reader, iplane1);
imwrite(img1, 'Colon1.tif');

%% step 2: train a classifier on the nuclei
% try to get the get nuclei completely but separe them where you can
% save as both simple segmentation and probabilities


%% step 3: use h5read to read your Ilastik simple segmentation
% and display the binary masks produced by Ilastik 
filename= 'Colon1_Simple Segmentation.h5';
datasetname = '/exported_data';
data = h5read(filename,datasetname);

% (datasetname = '/exported_data')
% Ilastik has the image transposed relative to matlab
% values are integers corresponding to segmentation classes you defined,
% figure out which value corresponds to nuclei
data_trans = transpose(squeeze(data)); 
data_trans(data_trans == 1) = 0;%1 is background  
data_trans(data_trans == 2 )= 1;%2 is nucleus


%% step 3.1: show segmentation as overlay on raw data
img_overlay = imfuse(data_trans, img1); 
imshow(img_overlay);

%imshow (cat(3, imadjust(img1), imadjust(data_trans), zeros(size(img1)))); 

%% step 4: visualize the connected components using label2rgb
% probably a lot of nuclei will be connected into large objects
rgb = label2rgb(bwlabel(data_trans), 'jet', 'w', 'shuffle'); 
imshow(rgb);
%% step 5: use h5read to read your Ilastik probabilities and visualize
filename1='Colon1_Probabilities.h5';
datasetname1 = '/exported_data';
data1 = h5read(filename1,datasetname1);
data_prob_chan1 = transpose(squeeze(data1(1,:,:))); %fish out channel 1 and squeeze
data_prob_chan2 = transpose(squeeze(data1(2,:,:))); %fish out channel 2 and squeeze
%no need for the following%
%data_prob_chan1(data_prob_chan1 == 1) = 0;  
%data_prob_chan1(data_prob_chan1 == 2 )= 1;
%data_prob_chan2(data_prob_chan2 == 1) = 0;  
%data_prob_chan2(data_prob_chan2 == 2 )= 1;

imshow(data_prob_chan1)
title('channel1')
figure
imshow(data_prob_chan2)
title('channel2')
% it will have a channel for each segmentation class you defined

%% step 6: threshold probabilities to separate nuclei better

%make a loop to see which one works the best
for i = 0.1:0.1:1.0 
    newdata = data_prob_chan2 < i;
    figure
    imshow(newdata);
end
% 0.5 gives a good result. 
probability = 0.5
newdata = data_prob_chan1 < probability;
imshow(newdata);
rgb_newdata = label2rgb(bwlabel(newdata), 'jet', 'w', 'shuffle'); 
imshow(rgb_newdata);
title('prob thresholding')
figure 
rgb = label2rgb(bwlabel(data_trans), 'jet', 'w', 'shuffle'); 
imshow(rgb);

%% step 7: watershed to fill in the original segmentation (~hysteresis threshold)
CC = bwconncomp(newdata);
stats = regionprops(CC, 'Area');
area = [stats.Area];
fusedCandidates = area>mean(area) + std(area); 
sublist = CC.PixelIdxList(fusedCandidates);
sublist = cat(1, sublist{:}); 
fusedMask = false(size(newdata));
fusedMask(sublist) = 1; 
imshow(fusedMask, 'InitialMagnification', 'fit')
title('fusedmask')
%eroding the centers
s = round(1.2*sqrt(mean(area))/pi);
nucmin = imerode(fusedMask,strel('disk',s));
imshow(nucmin,'InitialMagnification', 'fit');
title('nucmin') 
%getting the regions outside
outside = ~imdilate(fusedMask, strel('disk',1));
imshow(outside, 'InitialMagnification', 'fit');
title('outside') 
%define basin for watershed
basin = imcomplement(bwdist(outside));
basin= imimposemin(basin, nucmin | outside); 
L = watershed(basin); 
imshow(L); colormap('jet'); caxis([0 20]);
%combining mask
newMask = L>1 | (newdata - fusedMask); 
imshow(newMask, 'InitialMagnification', 'fit');
% my mask doesnot look impressive. 

%% step 8: perform hysteresis thresholding in Ilastik and compare the results
% explain the differences
% I cannot remember how to do this:( 

%% step 9: clean up the results more if you have time 
% using bwmorph, imopen, imclose etc
% Sorry I don't have time. 
