%% Jessie Talasco Final Project
%% 1.2 Visualization of Data
% Load the files and display all images in a sequential fashion
filename='Control.mat';
load(filename,'ImageData')
j=size(ImageData);
for i=1:j(3)
    imshow(ImageData(:,:,i))
    colorbar
end
fprintf('GEF.mat appears to be much more pixelated than the other two data files. I expect this will make the cell and membrane masks more difficult.')
%% 2.1 ROI Selection
% Crop the image to exclude the time stamp and make a figure showing the
% image before and after cropping at the first time point.
t=59;
cropped=ImageData(16:end,:,t);
figure
subplot(1,2,1), imshow(ImageData(:,:,1))
subplot(1,2,2), imshow(cropped)
%% 2.2 Thresholding to Select Membranes
% Apply thresholding technique to identify pixels belonging to membranes
% using 60 as the threshold value to make a membrane mask. Make a figure
% and discuss whether 60 is an appropriate value.
figure
threshold=60;
thresholded=cropped>threshold;
imshow(thresholded)
fprintf('60 seems to be an appropriate threshold there are clear white boundaries here that seem to make the membranes very distinctive.')
fprintf('If the threshold is drastically increased the membranes appear much larger than they actually are. There is lots of white space. If the threshold is increased dramatically the opposite would happen and you would miss some of the membrane.')
%% 2.3 Generation of Cell Mask Based on Membrane Mask
% Make a cell mask based on the membrane mask and make a figure
cellmask=~thresholded;
figure
subplot(1,2,1), imshow(thresholded)
colorbar
title(sprintf('Membranes, time=%dmin,threshold=%d',t-1,threshold))
subplot(1,2,2), imshow(cellmask)
colorbar
title(sprintf('Cells, time=%dmin,threshold=%d',t-1,threshold))
%% 2.4 Otsu's Method for Automatic Selection of Threshold
% Use Multithresh to find appropriate threshold, double check
% appropriateness by making a figure and discuss
thresh=multithresh(cropped);
thresholded2=cropped>thresh;
figure
subplot(1,2,1), imshow(thresholded2)
title(sprintf('Membranes, time=%dmin,threshold=%d',t-1,thresh))
colorbar
subplot(1,2,2), imshow(cropped)
title(sprintf('Original, time=%dmin,threshold=%d',t-1,threshold))
colorbar
%% 3.1 Quantification of Cell Area: Single Time Point
% Use bwconncomp and regionprops on the cell mask to generate figures close
% to figure 2(a), 2(c), and 2(e). Output mean and standard deviation of
% cell areas for a given image. Identify what went wrong in figure 2(b)
thresh=56;
thresholded2=cropped>thresh;
cellmask=~thresholded2;
CC=bwconncomp(cellmask,4);
stats=regionprops(CC,'Centroid','Area','PixelList');
figure
imshow(cellmask)
hold on
for k=1:length(stats)
    c=stats(k).Centroid;
    text(c(1), c(2), num2str(k), 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
figure;
imshow(cellmask)
hold on
[rows,cols]=size(cropped);
topedge=[1:cols; ones(1,cols)]';
bottomedge=[1:cols; rows * ones(1, cols)]';
left=[ones(1, rows); 1:rows]'; 
right=[cols * ones(1, rows); 1:rows]';
targets=[topedge; bottomedge; left; right];
n=0;
areas=[];
for k=1:length(stats)
    matches=ismember(stats(k).PixelList, targets, 'row');
    if~ any(matches)
         n=n+1;
         c=stats(k).Centroid; 
         text(c(1), c(2), num2str(n), 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
         areas=[areas ,stats(k).Area];
     end
end
figure;
h=histogram(areas);
xlim([0 max(h.BinEdges)])
xlabel('Area (# of pixels)')
ylabel('Frequency')
title('Histogram of Area: Control, time=0min')
m=mean(areas);
standd=std(areas);
fprintf('mean =%.2g\n,standard deviation =%.2g\n',m,standd)
fprintf('The threshold needed to be adjusted to better identifiy the cells using bwconncomp and regionprops. With too high of a threshold you end up with an image similar to figure 2(b) where multiple cells are seen as one and the centroids and areas are thrown off because the membrane boundary is not complete.')
%% 3.2 Temporal Changes of Cell Area Statistics
% make a function that performs a series ot tasks as in section 3.1 for any
% given 2D cell mask. Use the function to quantify and save the mean and
% standard deviation of cell areas for each time point for control data.
% Plot the coefficient of vatiation vs time point for control data.
j=size(ImageData);
ml=[];
standdl=[];
for i=1:j(3)
    cropped=ImageData(16:end,:,i);
    threshl=56;
    thresholdedl=cropped>threshl;
    cellmask=~thresholdedl;
    [m,standd]=CellAreaQuant(cellmask);
     ml=[ml ,m];
     standdl=[standdl, standd];
end
coeff=standdl./ml;
figure
plot(coeff,'-o', 'LineWidth', 2, 'MarkerSize', 6);
title('Coefficient of Variation vs Time')
xlabel('Time')
ylabel('Coefficient of Variation')
%% 4.0 Statistical Analysis
% Perform appropriate statistical test, find the p value, make a
% conclusion, create a figure to demonstrate the conclustion.
data=load('COV10min.txt');
data=data';
[p,tbl,stats]=anova1(data);
fprintf('P=%.2g\n, since this is less than alpha=0.05 we can conclude that not all of the means are the sane.',p)
