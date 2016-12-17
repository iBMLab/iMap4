function imwithpic(imagetmp,backgroundfile)

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
        im3D = repmat(rgb2gray(imbackground),[1,1,3]);
    else
        im3D = repmat(imbackground,[1,1,3]);
    end
end

im3D2 = imresize(im3D,[size(imagetmp,1),size(imagetmp,2)],'nearest');
range = [min(imagetmp(:)),max(imagetmp(:))];
toimagebgbeta = imagetmp.*0.7+im3D2.*0.3;
imshow(toimagebgbeta,range);hold on
% axis off

contv=linspace(min(imagetmp(:)),max(imagetmp(:)),6);
if isfinite(contv)
    imcontour(imagetmp,contv);colorbar;caxis([0 range(end)])
end
set(gca,'XTick',[],'YTick',[])


