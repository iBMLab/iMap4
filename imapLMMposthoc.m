function [Posthoc]=imapLMMposthoc(StatMap,FixMap,LMMmap,method,flag,formula2)
% Usage: [PostHoc]=imapLMMposthoc(StatMap,FixMap,LMMmap,method,flag,formula2)
%           method   - 'mean'/'sum' value in the cluster
%           flag     - 1 display result (default)
%           formula2 - using another LMM formula other than the original
%                      model to perform posthoc
% post-hoc contrast on raw/smoothed data (total fixation duration or fixation
% number), based on main effect cluster selected by hand.
% notice: mean fixation duration could be computed by total fixation
% duration./fixation number
% See also imapLMM, imapLMMcontrast, imapLMMdisplay
%
% 2015-02-12 Junpeng Lao, University of Fribourg.
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2015

scrsz=get(0,'ScreenSize');% get screen size for output display
tbl=LMMmap.Variables;
if isa(tbl,'dataset')
    VarNames=tbl.Properties.VarNames;
elseif isa(tbl,'table')
    VarNames=tbl.Properties.VariableNames;
    tbl= table2dataset(tbl);
end
if nargin<5
    flag=1;
end
if nargin<6
    formula=LMMmap.Formula;
else
    formula=formula2;
end
varargin=LMMmap.FitOptions;
lmetmp = LinearMixedModel.fit(tbl,char(formula),varargin{:});

DesignMatrix=designMatrix(lmetmp);

mask=StatMap.Pmask;
Label=StatMap.label;
mask2=zeros(size(mask));
ncluster=0;
nmask=size(mask,1);
for imask=1:nmask
    [lbmask,num]=bwlabel(squeeze(mask(imask,:,:)));
    mask2(imask,:,:)=lbmask+ncluster;
    ncluster=ncluster+num;
end
coefname=lmetmp.CoefficientNames;
categypredi=lmetmp.VariableInfo.InModel&lmetmp.VariableInfo.IsCategorical;
categyprediname=VarNames(categypredi);
% remove subject/grouping/random predictor if there is any
exclu1=zeros(length(categyprediname),1);
for icate=1:length(categyprediname)
    catename=categyprediname{icate};
    logictmp=strncmp(catename,coefname,length(catename));
    if sum(logictmp)==0
        exclu1(icate)=1;
    end
end
categyprediname(exclu1==1)=[];
Ncatepred=length(categyprediname);
if Ncatepred==0
    error('Categorical predictor not found, no post-hoc computation is performed')
else
    label=cell(length(tbl),Ncatepred);
    strlabel=cell(length(tbl),1);
    for icc=1:Ncatepred
        label(:,icc)=cellstr(eval(['tbl.' categyprediname{icc}]));
    end
    for ii=1:length(label)
        strlabel{ii,:}=strjoin(label(ii,:),'_');
    end
    [CatePredictor,btmp]=unique(strlabel,'rows');
    CateContrast=DesignMatrix(btmp,:);
    CateContrast(CateContrast(:)~=0&CateContrast(:)~=1&CateContrast(:)~=-1)=0;
    contrast=num2cell(CateContrast,2);
end

icinfo=1;
clusterinfo=[];
for imask=1:nmask
    h1=figure('NumberTitle','off','Name','Select the cluster or clusters for Post-hoc. Press Enter/Return to continue','Position',[1 1 scrsz(3)/2 scrsz(4)/2]);
    masktmp=squeeze(mask2(imask,:,:));
    imagesc(masktmp);hold on
    axis equal;axis off
    PositionCl=zeros(max(masktmp(:)),2);
    imm=1;
    xy=[-1 -1];
    while imm<(max(masktmp(:))+1) && isempty(xy)==0
        xy=ginput(1);
        if isempty(xy)==0
            PositionCl(imm,:) = xy;
            scatter(PositionCl(imm,1),PositionCl(imm,2),100,[1 1 1],'filled')
            imm=imm+1;
        end
    end
    close(h1)
    PositionCl(imm:end,:)=[];
    if isempty(PositionCl)==1
        disp('No cluster is selected from the current map')
    else
        PositionCl=round(PositionCl);
        for iposi=1:size(PositionCl,1)
            clusterinfo(icinfo,:)=[imask, masktmp(PositionCl(iposi,2),PositionCl(iposi,1))];
            icinfo=icinfo+1;
        end
    end
end

clusterinfo(clusterinfo(:,2)==0,:)=[];
if size(clusterinfo,1)>1
    disp('More than 1 cluster is selected, values within clusters will be combine')
end

maskfinal=zeros(size(masktmp));
for ii=1:size(clusterinfo,1)
    maskfinal(mask2(clusterinfo(ii,1),:,:)==clusterinfo(ii,2))=1;
end
maskfinal=logical(maskfinal);

switch method
    case 'mean'
        tbl.PixelIntensity=nanmean(FixMap(:,maskfinal),2);
    case 'sum'
        tbl.PixelIntensity=nansum(FixMap(:,maskfinal),2);
    otherwise
        error('Please input the method as mean or sum')
end

lmeposthoc = LinearMixedModel.fit(tbl,char(formula),varargin{:});

numcomp    = length(contrast);
posthocmat = NaN(numcomp,numcomp,3);
betao      = double(lmeposthoc.Coefficients(:,2));
betacov    = lmeposthoc.CoefficientCovariance;
alpha      = .05;
BetaCate   = array2table(zeros(numcomp,3),'VariableNames',{'beta' 'Lower' 'Upper'},'RowNames',CatePredictor);
dfre       = lmeposthoc.DFE;

for ia=1:numcomp
    BetaCate.beta(ia) = contrast{ia}*betao;
    delta = squeeze(tinv(1-alpha/2,dfre).*sqrt((contrast{ia}*betacov*contrast{ia}')))';
    BetaCate.Lower(ia)= BetaCate.beta(ia)-delta;
    BetaCate.Upper(ia)= BetaCate.beta(ia)+delta;
    for ib=(ia+1):numcomp
        contrst=contrast{ia}-contrast{ib};
        posthocmat(ia,ib,1)=contrst*betao;
        [posthocmat(ia,ib,3),posthocmat(ia,ib,2),df1,df2]=coefTest(lmeposthoc,contrst);%
    end
end
posthocmat(:,:,2) = sign(posthocmat(:,:,1)).*sqrt(posthocmat(:,:,2));
Posthoc.SelectCluster = maskfinal;
Posthoc.DF            = [df1,df2];
Posthoc.betadiff      = mat2dataset(posthocmat(:,:,1),'VarNames',CatePredictor,'ObsNames',CatePredictor);
Posthoc.Tval          = mat2dataset(posthocmat(:,:,2),'VarNames',CatePredictor,'ObsNames',CatePredictor);
Posthoc.pval          = mat2dataset(posthocmat(:,:,3),'VarNames',CatePredictor,'ObsNames',CatePredictor);
Posthoc.BetaCate      = BetaCate;
if flag==1
    %%
    figure('NumberTitle','off','Name','Post-hoc Mask','Position',[1 scrsz(4)/3 scrsz(3)/2 scrsz(4)/2]);
    imagesc(Posthoc.SelectCluster)
    axis equal;axis off
   
    figure('NumberTitle','off','Name','Post-hoc Statistics (Left minus Bottom)','Position',[1 1 scrsz(3) scrsz(4)/2]);
    warning('off')
    mat=double(Posthoc.betadiff);
    subplot(1,3,1)
    imsqrmat(mat, Posthoc.betadiff.Properties.VarNames,Posthoc.betadiff.Properties.ObsNames);
    title('Betas')
    
    mat=double(Posthoc.Tval);
    subplot(1,3,2)
    imsqrmat(mat, Posthoc.betadiff.Properties.VarNames,Posthoc.betadiff.Properties.ObsNames);
    title('Tvalue')
    
    mat2=double(Posthoc.pval);
    subplot(1,3,3)
    imagesc(mat2<(.05/sum(~isnan(mat2(:))))); 
    axis square;     
    set(gca,'XTick',1:size(mat2,1),...                         %# Change the axes tick marks 
        'XTickLabel',Posthoc.betadiff.Properties.VarNames,...  %#   and tick labels 
        'YTick',1:size(mat2,1),... 
        'YTickLabel',Posthoc.betadiff.Properties.ObsNames,... 
        'TickLength',[0 0],... 
        'xdir','reverse','ydir','normal'); 
    title('pValue < .05 (Bonferroni Corrected)')
end
