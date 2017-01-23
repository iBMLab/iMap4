function mapwithbg(image,backgroundfile,cmap,range,mask)
% heatmap on top of a background, with a 6 level contour map to highlight
% the value.
% input: image, backgroundfile
% additional:
%             cmap  - colormap, any input accepted by matlab colormap
%             function
%             range - map display range
%             mask  - display as a black contour line
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2016

if ischar(backgroundfile)&&~isempty(backgroundfile)
    imbackground = double(imread(sprintf(backgroundfile)))/255;
elseif isempty(backgroundfile)
    im3D=[];
    imbackground = [];
else
    imbackground = double(backgroundfile)/max(double(backgroundfile(:)));
end

if ~isempty(imbackground)
    % transfer background image to gray scale
    if size(imbackground,3)>1
        im3D2 = repmat(rgb2gray(imbackground),[1,1,3]);
    else
        im3D2 = repmat(imbackground,[1,1,3]);
    end
end

if nargin<3 || isempty(cmap)
    colormap('default');
    colourmap=colormap;
elseif ischar(cmap)&& ~isempty(cmap)
    cmap2=colormap(cmap);
    colourmap=cmap2;
end

if nargin<4 || isempty(range)
    range = [min(image(:)),max(image(:))];
end

if nargin<5
    mask = [];
end

ncontour = 6;

imagetmp      = imresize(image,[size(im3D2,1),size(im3D2,2)],'nearest');
toimagesg     = indtorgb(imagetmp,range(1),range(2),colourmap);
toimagebgbeta = toimagesg.*0.7+im3D2.*0.3;

imshow(toimagebgbeta,range);hold on
% axis off

contv=linspace(min(imagetmp(:)),max(imagetmp(:)),ncontour);
if isfinite(contv)
    imcontour(imagetmp,contv);colorbar;caxis(range)
end
if isfinite(mask)
    masktmp       = imresize(mask,[size(im3D2,1),size(im3D2,2)],'nearest');
    imcontour(1:size(masktmp,2),1:size(masktmp,1),masktmp,1,'k','LineWidth',1.5)                    
end
set(gca,'XTick',[],'YTick',[])


