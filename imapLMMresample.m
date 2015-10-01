function [ResampStat]=imapLMMresample(FixMap,LMMmap,c,h,effect,method,nboot,grouping,rmRE)
% imapLMMmontecarlo performs a nonparametric statistical test by calculating
% Monte-Carlo estimates of the significance probabilities and/or critical values
% from the resampling distribution.
% This function is called by imapLMMmcc with bootstrap or permutation
% option, but you can call it independently as well.
% input: FixMap, mask,
%       DX - fixed effect design matrix
%       Zx - random effect design matrix
%        c - contrast matrix
%        h - hypothesis matrix
%   effect - fixed/random
%   method - permutation/bootstrap
%    nboot - number of resampling
% grouping - specify group index to keeping the group variance constant
%     rmRE - 1 remove random effect, 0 keeping subject variance
% output: ResampStat with field {parameters} {resampleTABLE} {resampleFvalue}
%            {resamplePvalue} {resmapleBeta}
%
% 2014-09-14 Junpeng Lao, University of Fribourg.
% 2014-12-02 Option implemented for random with replacement within group
%            during bootstrap
%
% Copyright (C) iMap Team 2015
v=ver;
has_fsolve=any(strcmp({v.Name},'Parallel Computing Toolbox'));
if has_fsolve==1
    if isfield(LMMmap.runopt,'parallelname')==0
        gridname='local';
    else
        gridname=LMMmap.runopt.parallelname;
    end
end
warning('off')
effect=lower(effect);
Zx=LMMmap.RandomEffects.DX;
tbl=LMMmap.Variables;

% exclude empty trials/conditions
outtrial=full(sum(Zx,2))==0;
tbl(outtrial,:)=[];
Zx(outtrial,:)=[];
FixMap(outtrial,:,:)=[];
grouping(outtrial,:)=[];

if isa(tbl,'dataset')
    VarNames=tbl.Properties.VarNames;
elseif isa(tbl,'table')
    VarNames=tbl.Properties.VariableNames;
    tbl=table2dataset(tbl);
end
coefname=LMMmap.CoefficientNames;
categypredi=LMMmap.VariableInfo.InModel&LMMmap.VariableInfo.IsCategorical;
categyprediname=VarNames(categypredi);
contiupredi=LMMmap.VariableInfo.InModel&~LMMmap.VariableInfo.IsCategorical;
contiuprediname=VarNames(contiupredi);
% remove subject/grouping/random predictor if there is any
exclu1=zeros(length(categyprediname),1);
exclu2=zeros(length(contiuprediname),1);

for icate=1:length(categyprediname)
    if isempty(strmatch(categyprediname{icate},coefname))
        exclu1(icate)=1;
    end
end
categyprediname(exclu1==1)=[];
for iconu=1:length(contiuprediname)
    if isempty(strmatch(contiuprediname{iconu},coefname))
        exclu2(iconu)=1;
    end
end
contiuprediname(exclu2==1)=[];
label=cell(length(tbl),length(categyprediname));
if ~isempty(label)
    strlabel=cell(length(tbl),1);
    for icc=1:length(categyprediname)
        label(:,icc)=cellstr(eval(['tbl.' categyprediname{icc}]));
    end
    for ii=1:length(label)
        strlabel{ii,:}=strjoin(label(ii,:),'_');
    end
    [CatePredictor,~]=unique(strlabel,'rows');
end
[GroupPredictor,~]=unique(grouping,'rows');

if strcmp(effect,'random')
    effect2='random';
elseif strcmp(effect,'fixed') || strcmp(effect,'model beta')
    DX=LMMmap.modelX;
    effect2='fixed';
elseif strcmp(effect,'predictor beta')
    DX=LMMmap.SinglePred.DesignMatrix;
    effect2='fixed';
end
DFmodel=LMMmap.modelDFE;

mask=isnan(LMMmap.MSE)==0;
nonnan=find(mask==1);
[Nitem,cSize,rSize]=size(FixMap);
ResampStat=struct;
parameters{1}=c;
parameters{2}=h;
parameters{3}=effect;
parameters{4}=method;
parameters{5}=nboot;
ResampStat.params=parameters;

switch effect2
    case 'fixed'
        % prelocate memory
        Nc=length(c);
        
        DF=zeros(Nc,2);
        DF(:,2)=LMMmap.modelDFE;
        for ic=1:Nc
            DF(ic,1)=rank(c{ic});
        end
        
        if rmRE==1
            % for Y ~ X*b + Z*B + e, remove random effect Z*B
            ZB=NaN(Nitem,cSize,rSize);
            BLUP=squeeze(LMMmap.RandomEffects.RandomStat(:,1,:,:));
            for it=1:Nitem
                ZB(it,:)=Zx(it,:)*BLUP(:,:);
            end
            Y=FixMap-ZB;
        else
            Y=FixMap;
        end
        Yfixed=Y(:,mask);
        % Yfixed=FixMap(:,mask);
        
        % reproduce the beta (should be the same as fixed effect beta from
        % the LMM estimation)
        % beta=(reshape(DX\Y(:,:),size(DX,2),size(FixMap,2),size(FixMap,3)));
        
        % take data within mask.
        npixel=size(Yfixed,2);
        switch method
            case 'permutation'
                %%
                bY=Yfixed;
                Forg=NaN(Nc,cSize,rSize);
                resampleFvalue=NaN(nboot,Nc,cSize,rSize);
                for ic=1:Nc
                    c2=c{ic};
                    h2=h{ic};
                    
                    % Partition design matrix
                    cu=null(c2,'r');
                    % C=[c2;cu'];
                    c2=c2';
                    % the following part is from Appendix A. of Winkler, et al (2014).
                    % Permutation inference for the general linear model. Neuroimage, 92,
                    % 381-397.
                    D=pinv(DX'*DX);
                    cv=cu-c2*pinv(c2'*D*c2)*c2'*D*cu;
                    parX=DX*D*c2*pinv(c2'*D*c2);
                    parZ=DX*D*cv*pinv(cv'*D*cv);
                    M=[parX,parZ];
                    
                    cnew=zeros(size(c2'));
                    for ir=1:rank(parX)
                        cnew(ir,ir)=1;
                    end
                    
                    % contruct permutation table
                    % 2015-05-05
                    % Permutation table should be compute according to Winkler,
                    % et al (2014). Permutation inference for the general
                    % linear model. Neuroimage, 92, 381-397.
                    %
                    % We dont do sign flips
                    B=1;
                    boot_index1=zeros(nboot+500,Nitem);
                    while size(unique(boot_index1,'rows'),1)<nboot+1
                        tmp = randperm(Nitem);
                        boot_index1(B,:)=tmp;
                        B=B+1;
                    end
                    boot_index=unique(boot_index1(boot_index1(:,1)~=0,:),'rows');
                    ResampStat.resTABLE=boot_index;
                    
                    % orignial statistic for the fixed effect (revisit)
                    % bY2=(1-parZ*pinv(parZ))*bY;
                    obeta=M\bY;
                    cz=ones(1,rank(M))-sum(cnew,1);
                    bY2=bY-M(:,cz==1)*obeta(cz==1,:);

                    bbeta=M\bY2;
                    bMSE=sum((M*bbeta-bY2).^2)./DFmodel;
                    bMSEinv=bMSE.^-1;
                    % covprj=inv(M'*M);
                    % instead of calculating the covariance matrix and
                    % the inverse of the quadratic form c*covb*c'; here
                    % do QR decomposition to get the cholesky factor
                    X=qr(M,0);Rtmp=triu(X);
                    R=Rtmp(1:size(Rtmp,2),:);
                    S=inv(R);
                    
                    quadformCOV=cnew/R*S'*cnew';%same as (c*(S*S')*c) == (c*covb*c')
                    qualCOVinv=inv(quadformCOV);
                    cbeta_h=(cnew*bbeta-h2)';
                    
                    I=repmat(1:npixel,[DF(ic,1),1]);
                    trc=cbeta_h';
                    cmt_sp = sparse(I(:),1:npixel*DF(ic,1),trc);
                    
                    covb=repmat(qualCOVinv,[1,1,npixel]);
                    I = repmat(reshape(1:DF(ic,1)*npixel,DF(ic,1),1,npixel),[1 DF(ic,1) 1]);
                    J = repmat(reshape(1:DF(ic,1)*npixel,1,DF(ic,1),npixel),[DF(ic,1) 1 1]);
                    covb_sp = sparse(I(:),J(:),covb(:));
                    
                    tmp_sp=cmt_sp*covb_sp;
                    tmpproj=tmp_sp*cmt_sp';
                    Itmp=speye(npixel);
                    proj=tmpproj(logical(Itmp));
                    
                    Forg(ic,nonnan)=bMSEinv'.*proj./DF(ic,1);
                    
                    % resampling nboot times
                    if has_fsolve==1
                        % parpool;
                        try 
                            parpool(gridname); 
                            pctRunOnAll warning('off')
                        end
                        parfor ib=1:nboot
                            % contruct index (we do this in this loop to adapt for unbalance design)
                            bs=boot_index(ib,:);
                            bDX=M(bs,:);
                            
                            % orignial statistic for the fixed effect (revisit)
                            bbeta=bDX\bY2;
                            bMSE=sum((bDX*bbeta-bY2).^2)./DFmodel;
                            bMSEinv=bMSE.^-1;
                            % covprj=inv(bDX'*bDX);
                            % instead of calculating the covariance matrix and
                            % the inverse of the quadratic form c*covb*c'; here
                            % do QR decomposition to get the cholesky factor
                            X=qr(bDX,0);Rtmp=triu(X);
                            R=Rtmp(1:size(Rtmp,2),:);
                            S=inv(R);
                            
                            quadformCOV=cnew/R*S'*cnew';%same as (c*(S*S')*c) == (c*covb*c')
                            qualCOVinv=inv(quadformCOV);
                            cbeta_h=(cnew*bbeta-h2)';
                            
                            I=repmat(1:npixel,[DF(ic,1),1]);
                            trc=cbeta_h';
                            cmt_sp = sparse(I(:),1:npixel*DF(ic,1),trc);
                            
                            covb=repmat(qualCOVinv,[1,1,npixel]);
                            I = repmat(reshape(1:DF(ic,1)*npixel,DF(ic,1),1,npixel),[1 DF(ic,1) 1]);
                            J = repmat(reshape(1:DF(ic,1)*npixel,1,DF(ic,1),npixel),[DF(ic,1) 1 1]);
                            covb_sp = sparse(I(:),J(:),covb(:));
                            
                            tmp_sp=cmt_sp*covb_sp;
                            tmpproj=tmp_sp*cmt_sp';
                            Itmp=speye(npixel);
                            proj=tmpproj(logical(Itmp));
                            
                            Ftmp1=bMSEinv'.*proj./DF(ic,1);
                            resampleFvalue(ib,ic,nonnan)=Ftmp1;
                        end
                        % delete(gcp)
                    else
                        waith=waitbar(0,'Resampling...');
                        for ib=1:nboot
                            waitbar(ib / nboot)
                            % contruct index (we do this in this loop to adapt for unbalance design)
                            bs=boot_index(ib,:);
                            bDX=M(bs,:);
                            
                            % orignial statistic for the fixed effect (revisit)
                            bbeta=bDX\bY2;
                            bMSE=sum((bDX*bbeta-bY2).^2)./DFmodel;
                            bMSEinv=bMSE.^-1;
                            % covprj=inv(bDX'*bDX);
                            % instead of calculating the covariance matrix and
                            % the inverse of the quadratic form c*covb*c'; here
                            % do QR decomposition to get the cholesky factor
                            X=qr(bDX,0);Rtmp=triu(X);
                            R=Rtmp(1:size(Rtmp,2),:);
                            S=inv(R);
                            
                            quadformCOV=cnew/R*S'*cnew';%same as (c*(S*S')*c) == (c*covb*c')
                            qualCOVinv=inv(quadformCOV);
                            cbeta_h=(cnew*bbeta-h2)';
                            
                            I=repmat(1:npixel,[DF(ic,1),1]);
                            trc=cbeta_h';
                            cmt_sp = sparse(I(:),1:npixel*DF(ic,1),trc);
                            
                            covb=repmat(qualCOVinv,[1,1,npixel]);
                            I = repmat(reshape(1:DF(ic,1)*npixel,DF(ic,1),1,npixel),[1 DF(ic,1) 1]);
                            J = repmat(reshape(1:DF(ic,1)*npixel,1,DF(ic,1),npixel),[DF(ic,1) 1 1]);
                            covb_sp = sparse(I(:),J(:),covb(:));
                            
                            tmp_sp=cmt_sp*covb_sp;
                            tmpproj=tmp_sp*cmt_sp';
                            Itmp=speye(npixel);
                            proj=tmpproj(logical(Itmp));
                            
                            Ftmp1=bMSEinv'.*proj./DF(ic,1);
                            resampleFvalue(ib,ic,nonnan)=Ftmp1;
                        end
                        close(waith)
                    end
                end
                ResampStat.Forg=Forg;
                ResampStat.resFvalue=resampleFvalue;
            case 'bootstrap'
                %%
                resampleFvalue=NaN(nboot,Nc,cSize,rSize);
                resamplePvalue=NaN(nboot,Nc,cSize,rSize);
                resmapleBeta=NaN(nboot,Nc,cSize,rSize);
                % center data
                % beta=DX\Yfixed;
                % Yc=Yfixed-DX*beta;
                Yc=zeros(size(Yfixed));
                if ~isempty(label)
                    for icate=1:length(CatePredictor)
                        tmpcate=CatePredictor{icate};
                        c_ind=strcmp(strlabel,tmpcate);
                        Yc(c_ind,:)=Yfixed(c_ind,:)-repmat(mean(Yfixed(c_ind,:)),[sum(c_ind),1]);
                    end
                else
                    Yc=Yfixed(randperm(Nitem,Nitem),:);
                end
                Yc=Yc(randperm(Nitem,Nitem),:);
                
                % contruct bootstrap table
                % random with replacement for all subject
                % find subject vector
                [~,b]=find(Zx~=1&Zx~=0);
                if ~isempty(b)
                    Zx(:,b)=[];
                end
                
                indxsbj=full(Zx);
                if sum(indxsbj(:)==1)==Nitem
                    indxsbj(Zx==1)=1:Nitem;
                else
                    try
                        Nmatrix=round(sum(indxsbj(:)==1)/Nitem);
                        columnsum=sum(indxsbj,1);
                        matrixblock=zeros(size(columnsum));
                        tmpsum=0;blocktype=1;
                        for ii=1:size(columnsum,2)
                            tmpsum=tmpsum+columnsum(ii);
                            if tmpsum<=Nitem
                                matrixblock(ii)=blocktype;
                            elseif tmpsum>Nitem
                                blocktype=blocktype+1;
                                tmpsum=0;
                                matrixblock(ii)=blocktype;
                            end
                        end
                        NNmatrix=zeros(Nmatrix,1);
                        for iN=1:Nmatrix
                            NNmatrix(iN)=sum(matrixblock==iN);
                        end
                        sbjindxtmp=matrixblock==find(NNmatrix==max(NNmatrix));
                        indxsbj(:,~sbjindxtmp)=[];
                        indxsbj(indxsbj==1)=1:Nitem;
                    catch
                        error('imapLMMresample error: Can not retrieve Subject column for resampling.')
                    end
                end
                Ns=rank(indxsbj); % number of subject
                
                % sample with replacement, create boot_table
                % boot_index=randi(Ns,nboot,Ns);
                B=1;
                boot_index=zeros(nboot,Ns);
                while B~=nboot+1
                    if length(GroupPredictor)>1 % more than 2 groups
                        tmp = [];
                        for ig=1:length(GroupPredictor)
                            sbjtmp=find(sum(indxsbj(grouping==GroupPredictor(ig),:))~=0);
                            tmp=[tmp sbjtmp(randi(length(sbjtmp),1,length(sbjtmp)))];
                        end
                        boot_index(B,:) = tmp;
                        B=B+1;
                    else
                        tmp = randi(Ns,1,Ns);
                        if length(unique(tmp))>= 2 % at least 3 different observations per boot
                            boot_index(B,:) = tmp;
                            B=B+1;
                        end
                    end
                end
                ResampStat.resTABLE=boot_index;
                
                % resampling nboot times
                if has_fsolve==1
                    % parpool;
                    try
                        parpool(gridname);
                        pctRunOnAll warning('off')
                    end
                    parfor ib=1:nboot
                        % contruct index (we do this in this loop to adapt for unbalance design)
                        bindx=indxsbj(:,boot_index(ib,:));
                        bs=bindx(bindx~=0);
                        bDX=DX(bs,:);
                        bY=Yc(bs,:);
                        df2=size(DX,1)-size(DX,2);
                        
                        bbeta=bDX\bY;
                        bMSE=sum((bDX*bbeta-bY).^2)./DFmodel;
                        bMSEinv=bMSE.^-1;
                        % covprj=inv(bDX'*bDX);
                        % instead of calculating the covariance matrix and
                        % the inverse of the quadratic form c*covb*c'; here
                        % do QR decomposition to get the cholesky factor
                        X=qr(bDX,0);Rtmp=triu(X);
                        R=Rtmp(1:size(Rtmp,2),:);
                        S=inv(R);
                        
                        Ftmp1=NaN(Nc,npixel);
                        ptmp1=NaN(Nc,npixel);
                        betamap=NaN(Nc,npixel);
                        for ic=1:Nc
                            c2=c{ic};
                            % h2=h{ic};
                            h2=0;%%%% important as c*beta is already equal to 0 as the H0
                            quadformCOV=c2/R*S'*c2';%same as (c*(S*S')*c) == (c*covb*c')
                            qualCOVinv=inv(quadformCOV);
                            cbeta_h=(c2*bbeta-h2)';
                            if size(c2,1)>1
                                betamap(ic,:)=mean(cbeta_h,2);
                            else
                                betamap(ic,:)=c2*bbeta;
                            end
                            
                            I=repmat(1:npixel,[DF(ic,1),1]);
                            trc=cbeta_h';
                            cmt_sp = sparse(I(:),1:npixel*DF(ic,1),trc);
                            
                            covb=repmat(qualCOVinv,[1,1,npixel]);
                            I = repmat(reshape(1:DF(ic,1)*npixel,DF(ic,1),1,npixel),[1 DF(ic,1) 1]);
                            J = repmat(reshape(1:DF(ic,1)*npixel,1,DF(ic,1),npixel),[DF(ic,1) 1 1]);
                            covb_sp = sparse(I(:),J(:),covb(:));
                            
                            tmp_sp=cmt_sp*covb_sp;
                            tmpproj=tmp_sp*cmt_sp';
                            Itmp=speye(npixel);
                            proj=tmpproj(logical(Itmp));
                            
                            Ftmp1(ic,:)=bMSEinv'.*proj./DF(ic,1);
                            ptmp1(ic,:)=1-fcdf(Ftmp1(ic,:),DF(ic,1),df2);
                        end
                        resampleFvalue(ib,:,nonnan)=Ftmp1;
                        resamplePvalue(ib,:,nonnan)=ptmp1;
                        resmapleBeta(ib,:,nonnan)=betamap;
                    end
                    % delete(gcp)
                else
                    waith=waitbar(0,'Resampling...');
                    for ib=1:nboot
                        waitbar(ib / nboot)
                        % contruct index (we do this in this loop to adapt for unbalance design)
                        bindx=indxsbj(:,boot_index(ib,:));
                        bs=bindx(bindx~=0);
                        bDX=DX(bs,:);
                        bY=Yc(bs,:);
                        df2=size(DX,1)-size(DX,2);
                        
                        bbeta=bDX\bY;
                        bMSE=sum((bDX*bbeta-bY).^2)./DFmodel;
                        bMSEinv=bMSE.^-1;
                        % covprj=inv(bDX'*bDX);
                        % instead of calculating the covariance matrix and
                        % the inverse of the quadratic form c*covb*c'; here
                        % do QR decomposition to get the cholesky factor
                        X=qr(bDX,0);Rtmp=triu(X);
                        R=Rtmp(1:size(Rtmp,2),:);
                        S=inv(R);
                        
                        Ftmp1=NaN(Nc,npixel);
                        ptmp1=NaN(Nc,npixel);
                        betamap=NaN(Nc,npixel);
                        for ic=1:Nc
                            c2=c{ic};
                            % h2=h{ic};
                            h2=0;%%%% important as c*beta is already equal to 0 as the H0
                            quadformCOV=c2/R*S'*c2';%same as (c*(S*S')*c) == (c*covb*c')
                            qualCOVinv=inv(quadformCOV);
                            cbeta_h=(c2*bbeta-h2)';
                            if size(c2,1)>1
                                betamap(ic,:)=mean(cbeta_h,2);
                            else
                                betamap(ic,:)=c2*bbeta;
                            end
                            
                            I=repmat(1:npixel,[DF(ic,1),1]);
                            trc=cbeta_h';
                            cmt_sp = sparse(I(:),1:npixel*DF(ic,1),trc);
                            
                            covb=repmat(qualCOVinv,[1,1,npixel]);
                            I = repmat(reshape(1:DF(ic,1)*npixel,DF(ic,1),1,npixel),[1 DF(ic,1) 1]);
                            J = repmat(reshape(1:DF(ic,1)*npixel,1,DF(ic,1),npixel),[DF(ic,1) 1 1]);
                            covb_sp = sparse(I(:),J(:),covb(:));
                            
                            tmp_sp=cmt_sp*covb_sp;
                            tmpproj=tmp_sp*cmt_sp';
                            Itmp=speye(npixel);
                            proj=tmpproj(logical(Itmp));
                            
                            Ftmp1(ic,:)=bMSEinv'.*proj./DF(ic,1);
                            ptmp1(ic,:)=1-fcdf(Ftmp1(ic,:),DF(ic,1),df2);
                        end
                        resampleFvalue(ib,:,mask)=Ftmp1;
                        resamplePvalue(ib,:,mask)=ptmp1;
                        resmapleBeta(ib,:,mask)=betamap;
                    end
                    close(waith)
                end
                ResampStat.resFvalue=resampleFvalue;
                ResampStat.resPvalue=resamplePvalue;
                ResampStat.resBeta=resmapleBeta;
        end
        
    case 'random'
        error('Resampling method for random effect is not available yet')
end
warning('on')
end