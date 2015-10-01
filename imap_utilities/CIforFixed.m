function masknew=CIforFixed(StatMap,baseline)
% Confidence Interval for fixed effect beta
% size(StatMap_c.betaCI)
% size(StatMap_c.resampMat.resBeta)
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015
cmax=max(StatMap.betaCI(:));
cmin=min(StatMap.betaCI(:));
scrsz=get(0,'ScreenSize');% get screen size for output display

alpha=.05;
if isfield(StatMap,'resampMat')
    [nboot,Nc,xSize,ySize]=size(StatMap.resampMat.resBeta);
    masknew=zeros(2,Nc,xSize,ySize);
else
    [Nc,xSize,ySize]=size(StatMap.beta);
    masknew=zeros(1,Nc,xSize,ySize);
end
for ic=1:Nc
    betaORI=squeeze(StatMap.beta(ic,:,:));
    figure('Numbertitle','off','Name',StatMap.label{ic},'Position',[1 1 scrsz(3) scrsz(4)/2]);
    betaCI1=squeeze(StatMap.betaCI(ic,1,:,:));
    betaCI2=squeeze(StatMap.betaCI(ic,2,:,:));
    subplot(2,3,4)
    imagesc(betaCI1,[cmin cmax]);axis('equal','off')
    title('lower CI')
    subplot(2,3,1)
    imagesc(betaCI2,[cmin cmax]);axis('equal','off')
    title('upper CI')
    if isfield(StatMap,'resampMat')
        % compute Bootstrap CI
        betatmp=squeeze(StatMap.resampMat.resBeta(:,ic,:,:));
        betasort=sort(betatmp(:,:),1);
        % betaCI1b=reshape(betasort(round(nboot*alpha),:),[xSize,ySize])+betaORI;
        % betaCI2b=reshape(betasort(round(nboot*(1-alpha)),:),[xSize,ySize])+betaORI;
        betaCI1b=reshape(betasort(1,:),[xSize,ySize])+betaORI;
        betaCI2b=reshape(betasort(end,:),[xSize,ySize])+betaORI;
        subplot(2,3,5)
        imagesc(betaCI1b,[cmin cmax]);axis('equal','off')
        title('lower CI (bs)')
        subplot(2,3,2)
        imagesc(betaCI2b,[cmin cmax]);axis('equal','off')
        title('upper CI (bs)')
    end
    subplot(2,3,3)
    maskPM=betaCI1>baseline;
    imagesc(betaORI.*maskPM);colorbar;axis('equal','off')
    title('above average activation (Parametric CI)')
    % imagesc(betaCI1b-betaCI1);colorbar;axis('equal','off')
    % title('diff CI (bs-pm)')
    masknew(1,ic,:,:)=maskPM;
    if isfield(StatMap,'resampMat')
        subplot(2,3,6)
        maskBT=betaCI1b>baseline;
        imagesc(betaORI.*maskBT);colorbar;axis('equal','off')
        title('above average activation (Bootstrap CI)')
        masknew(2,ic,:,:)=maskBT;
    end
    % imagesc(betaCI2b-betaCI2);colorbar;axis('equal','off')
    % title('diff CI (bs-pm)')
end
