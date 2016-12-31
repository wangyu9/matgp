function [] = image_white2none(filenamein,filenameout,varargin)

% example useage: image_white2none('2.png','2t.png');

% read in .png file with alpha layer

method1 = false;

threshold = 255;
if(numel(varargin)>0)
   threshold = varargin{1}; 
   if(numel(varargin)>1)
        if(varargin{2}==true)
            method1 = true;
        end
   end
end

[img,~,alpha] = imread(filenamein);
  
alpha = ones(size(img,1),size(img,2));

if(method1)
    % not work do not know why rgb_sum = double(img(:,:,1))+double(img(:,:,2))+double(img(:,:,3));
    % alpha( find(rgb_sum)>=threshold ) = 0;

    alpha( find( (img(:,:,1)>=threshold)&(img(:,:,2)>=threshold)&(img(:,:,3)>=threshold) ) ) = 0;
    %alpha( find( img(:,:,1)==255&img(:,:,2)==255&img(:,:,3)==255 ) ) = 0;
else
    ksize = 10;
    
    A = threshold*ones(size(img,1),size(img,2));
    B = ones(ksize,ksize) ./ (ksize*ksize);
    
    S = conv2(A,B,'same');
    
    img2 = img;
    img2(:,:,1) = conv2(double(255-img(:,:,1)),B,'same');% because matlab always padded 0.
    img2(:,:,2) = conv2(double(255-img(:,:,2)),B,'same');
    img2(:,:,3) = conv2(double(255-img(:,:,3)),B,'same');
    
    alpha( find( 255-img2(:,:,1)>=threshold&255-img2(:,:,2)>=threshold&255-img2(:,:,3)>=threshold ) ) = 0;
end



imwrite(img,filenameout,'Alpha',alpha);