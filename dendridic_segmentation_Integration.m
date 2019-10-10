%% Pull file
% ref: http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/
Path = 'Images/Helicobacter_Pylori_Positive/RA035_1/';
File = 'RA035_1_Composite with Selection';
FileTif = strcat(Path,File,'.tif');
InfoImage = imfinfo(FileTif);
mImage = InfoImage(1).Width;
nImage = InfoImage(1).Height;
NumberImages = length(InfoImage);

%storage for image layers
FinalImage1 = zeros(nImage, mImage, NumberImages, 'uint16');

tifLink = Tiff(FileTif, 'r');
for i = 1:NumberImages
    tifLink.setDirectory(i);
    FinalImage1(:,:,i) = tifLink.read();
end
tifLink.close();

% this code is used in case multiple tiff files are used instead of one
% 
% File = '#34 a7-1 20x-4 thresholded nuclei';
% FileTif = strcat(Path,File,'.tif');
% InfoImage = imfinfo(FileTif);
% mImage = InfoImage(1).Width;
% nImage = InfoImage(1).Height;
% NumberImages = length(InfoImage);
% 
% FinalImage2 = zeros(nImage, mImage, NumberImages, 'uint16');
% 
% tifLink = Tiff(FileTif, 'r');
% for i = 1:NumberImages
%     tifLink.setDirectory(i);
%     FinalImage2(:,:,i) = tifLink.read();
% end
% tifLink.close();
% 
% %just to adjust for two .tifs instead of one with all three images
% [m,n] = size(FinalImage1);
% FinalImage = zeros(m,n,2);
% FinalImage(:,:,1) = FinalImage1;
% FinalImage(:,:,2) = FinalImage2;

FinalImage = FinalImage1;

%% Pulling Selection

%these need to be set manually --- changed based on image, would be a good
%place for user selection
% Nindex - layer for Nuclei Information
% DCindex - layer for Dendridic Information
Nindex = 1;
DCindex = 2;

% requires ROI to be saved as .csv from Fiji
SelectionPath = strcat(Path,"Selection.csv");
Selection = selection_logical(SelectionPath);

%padd ones if size doesn't match and error is thrown
Selection = padarray(Selection',[abs(size(FinalImage(:,:,Nindex),2)-size(Selection,2)) 2],1,'post')';
Selection = padarray(Selection',[abs(size(FinalImage(:,:,Nindex),1)-size(Selection,1)) 1],1,'post')';

% adjust size if ROI doesn't match image size
 Selection = Selection(1:size(FinalImage(:,:,Nindex),1),1:size(FinalImage(:,:,Nindex),2));



%% Nuclei Centers

% Segmentation of nuclei using DoG filter
gaussian1 = fspecial('Gaussian',[25,1],5);
gaussian2 = fspecial('Gaussian',[25,1],6);
Nuclei = convn(convn((FinalImage(:,:,Nindex))',gaussian1,'same')',gaussian1,'same') - convn(convn((FinalImage(:,:,Nindex))',gaussian2,'same')',gaussian2,'same') ;

%   generate a 'in process' image
% figure; imagesc(Nuclei(1:400,700:1100)); axis off;

% Finding Centers of Nuclei

Nuclei_Centers = imregionalmax(Nuclei);

%   generate a 'in process' image
% se = strel('disk',3);
% NC_Display = imdilate(Nuclei_Centers,se);
% figure; imagesc(NC_Display(1:400,700:1100)); axis off;

% pull off items outside selection area
Nuclei_Centers(~logical(full(Selection))) = 0;

%   generate a 'in process' image
% se = strel('disk',3);
% NC_Display = imdilate(Nuclei_Centers,se);
% figure; imagesc(NC_Display(1:400,700:1100)); axis off;

% correction for imregionalmax
Nuclei_Centers(Nuclei < 0.05*median(median(FinalImage(:,:,Nindex))) | FinalImage(:,:,Nindex) < 0.5*median(median(FinalImage(:,:,Nindex)))) = 0;

%   Results images
se = strel('disk',3);
NC_Display = imdilate(Nuclei_Centers,se);
%figure; imagesc(NC_Display); axis off;
     
    figure; imagesc(FinalImage(:,:,Nindex)); title("Nuclei Image"); 
    figure; imagesc(Nuclei); title("Nuclei after Convolution"); 
    figure; imagesc(NC_Display); title("Nuclei Center locations");
    figure; imagesc(150*NC_Display + 0.1*double(FinalImage(:,:,Nindex))); title("Nuclei over Nuclei Image");
    figure; imshowpair(FinalImage(:,:,Nindex),NC_Display); 
    
%% Filter dendridic cells

 gaussian3 = fspecial('Gaussian',15,2);
 
 Dendridic = convn(FinalImage(:,:,DCindex),gaussian3,'same');

% Threshold nuclei centers based on Dendridic Cells

 Nuclei_Dendridic = Nuclei_Centers;
 
 for k = [-2, -1, 0, 1, 2 ]

     tol = median(median(Dendridic)) + k*std(std(Dendridic));

     Nuclei_Dendridic(Dendridic(Nuclei_Centers == 0)<tol) = 0;


    % Display Center Locations
    xStart = 1;
    xEnd = 750;
    yStart = 1;
    yEnd = 750;
    se = strel('disk',3);
    ND_Display = imdilate(Nuclei_Dendridic,se);
    NC_Display = imdilate(Nuclei_Centers,se);
    figure; 
    subplot(2,4,1); imagesc(FinalImage(xStart:xEnd,yStart:yEnd,Nindex)); title("Nuclei Image"); 
    subplot(2,4,2); imagesc(Nuclei(xStart:xEnd,yStart:yEnd)); title("Nuclei after Convolution"); 
    subplot(2,4,3); imagesc(NC_Display(xStart:xEnd,yStart:yEnd)); title("Nuclei Center locations");
    subplot(2,4,4); imagesc(150*NC_Display(xStart:xEnd,yStart:yEnd) + 0.1*double(FinalImage(xStart:xEnd,yStart:yEnd,Nindex))); title("Nuclei over Nuclei Image");
    subplot(2,4,5); imagesc(FinalImage(xStart:xEnd,yStart:yEnd,DCindex)); title("Dendridic (D) Image"); 
    ax = gca; 
    clims = ax.CLim;
    subplot(2,4,6); imagesc(Dendridic(xStart:xEnd,yStart:yEnd)); title("D after Convolution"); 
    caxis(clims);
    subplot(2,4,7); imagesc(ND_Display(xStart:xEnd,yStart:yEnd)); title(strcat('median + ',num2str(k),' *std'));
    subplot(2,4,8); imagesc(150*ND_Display(xStart:xEnd,yStart:yEnd) + 0.1*double(FinalImage(xStart:xEnd,yStart:yEnd,DCindex))); title("D Nuclei over D Image");
    
    %save variables for Stats processing
    save(strcat('34c9_2_med_',num2str(k),'std'),'Selection','Nuclei_Centers','Nuclei_Dendridic');
 end

 


