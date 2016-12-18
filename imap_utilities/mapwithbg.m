function mapwithbg(image,backgroundfile,cmap,range,mask)
% map with background

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

if nargin<4
    range = [min(image(:)),max(image(:))];
end

if nargin<5
    mask = [];
end

imagetmp      = imresize(image,[size(im3D2,1),size(im3D2,2)],'nearest');
toimagesg     = indtorgb(imagetmp,range(1),range(2),colourmap);
toimagebgbeta = toimagesg.*0.7+im3D2.*0.3;

imshow(toimagebgbeta,range);hold on
% axis off

contv=linspace(min(imagetmp(:)),max(imagetmp(:)),6);
if isfinite(contv)
    imcontour(imagetmp,contv);colorbar;caxis([0 range(end)])
end
if isfinite(mask)
    imcontour(1:size(mask,2),1:size(mask,1),mask,1,'k','LineWidth',1)                    
end
set(gca,'XTick',[],'YTick',[])


