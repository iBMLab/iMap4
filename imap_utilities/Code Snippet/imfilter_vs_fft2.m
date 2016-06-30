%--------------------------------------------------------------------------
% code snippet to demonstrate smoothing (comparing speed).
% Copyright (C) iMap Team 2015
%%
% Generate a reasonable size image
% Notice: image larger than 500 pixel imfilter becomes too slow to apply.
xSize  = 400;
ySize  = 400;
Nfix   = xSize/10;
coordX = randi(xSize,Nfix,1);
coordY = randi(xSize,Nfix,1);

rawmap = full(sparse(coordY,coordX,1,ySize,xSize));

smoothingpic = 10;
[x, y] = meshgrid(-floor(xSize/2)+.5:floor(xSize/2)-.5, ...
                 -floor(ySize/2)+.5:floor(ySize/2)-.5);
gaussienne = exp(- (x .^2 / smoothingpic ^2) - (y .^2 / smoothingpic ^2));
gaussienne = (gaussienne - min(gaussienne(:))) ...
           / (max(gaussienne(:)) - min(gaussienne(:)));
f_fil      = fft2(gaussienne);

% fft2
f_mat = fft2(rawmap); % 2D fourrier transform on the points matrix
filtered_mat = f_mat .* f_fil;
smoothpic1 = real(fftshift(ifft2(filtered_mat))); % take the real part of the complex values from the fourier transform
f1 = @() real(fftshift(ifft2(fft2(rawmap) .* f_fil)));
fftime=timeit(f1);
% imfilter
smoothpic2 = imfilter(rawmap,gaussienne,'replicate','same','conv');
f2 = @() imfilter(rawmap,gaussienne,'replicate','same','conv');
imfiltime=timeit(f2);
% conv2
smoothpic3 = conv2(rawmap, gaussienne,'same');
f3 = @() conv2(rawmap, gaussienne,'same');
conv2time=timeit(f3);

% 
figure('Numbertitle','off');
subplot(2,3,1)
imagesc(gaussienne);axis square off
title('Gaussian Filter')
subplot(2,3,2)
scatter(coordX,coordY); set(gca,'YDir','reverse');axis square
title('Raw data - top left as (0,0)')
subplot(2,3,3);hold on
text(0,3,['fft2: ',num2str(fftime),'s']        ,'HorizontalAlignment','right','Interpreter','latex');
text(0,2,['imfilter: ',num2str(imfiltime),'s'] ,'HorizontalAlignment','right','Interpreter','latex');
text(0,1,['conv2 :',num2str(conv2time),'s']    ,'HorizontalAlignment','right','Interpreter','latex');
axis([-5,0,0,5],'off')
subplot(2,3,4)
imagesc(smoothpic1);axis square off
title('fft2 output')
subplot(2,3,5)
imagesc(smoothpic2);axis square off
title('imfilter output')
subplot(2,3,6)
imagesc(smoothpic3);axis square off
title('conv2 output')
