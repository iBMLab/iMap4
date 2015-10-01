%% descriptive result
load('DescriptvM_single_trial.mat')
load('FixMap_single_trial_scaled.mat')
%% filtering
% indx=find(DescriptvM.FixNum>200);
% FixMap(indx,:,:)=[];
% RawMap(indx,:,:)=[];
% 
% DescriptvM(indx,:)=[];
% PredictorM(indx,:)=[];
%%
ConditionM=DescriptvM(:,[6:end]);
MeasureM=DescriptvM(:,[1:5]);
nbins=50;
% overall distribution
figure('NumberTitle','off','Name','Eye movement measurement distribution');
% given that fixation number could be consider as a poisson process, the
% parameters could be fit with a Gamma distribution. Fixation
% duration/path length could both be consider some how as the waiting time
% of the possion process, thus should be appropriate to fit also a gamma
% distribution
subplot(3,2,1)
pd1=histfitimap(MeasureM.FixNum,nbins,'gamma');
title('Fixation Number')

subplot(3,2,3)
pd2=histfitimap(MeasureM.sumFixDur,nbins,'gamma');
title('Fixation Duration(SUM)')

subplot(3,2,4)
pd3=histfitimap(MeasureM.meanFixDur,nbins,'gamma');
title('Fixation Duration(MEAN)')

subplot(3,2,5)
pd4=histfitimap(MeasureM.totalPathLength,nbins,'gamma');
title('Path Length(TOTAL)')

subplot(3,2,6)
pd5=histfitimap(MeasureM.meanPathLength,nbins,'gamma');
title('Path Length(MEAN)')

%% box plot for each catigorical condition
CName=ConditionM.Properties.VarNames;
for ic=1:length(CName)
    indextmp=eval(['ConditionM.' CName{ic}]);
    figure('NumberTitle','off','Name',['Condition ' CName{ic}]);
    subplot(3,2,1)
    boxplot(MeasureM.FixNum,indextmp)
    title('Fixation Number')
    
    subplot(3,2,3)
    boxplot(MeasureM.sumFixDur,indextmp)
    title('Fixation Duration(SUM)')
    
    subplot(3,2,4)
    boxplot(MeasureM.meanFixDur,indextmp)
    title('Fixation Duration(MEAN)')
    
    subplot(3,2,5)
    boxplot(MeasureM.totalPathLength,indextmp)
    title('Path Length(TOTAL)')
    
    subplot(3,2,6)
    boxplot(MeasureM.meanPathLength,indextmp)
    title('Path Length(MEAN)')
end
%% mean fixation map
maxsub=20;
for ic=1:length(CName)
    %%
    h(1)=figure('NumberTitle','off','Name',['Mean fixation bias - ' CName{ic}]);
    conditiontmp=eval(['ConditionM.' CName{ic}]);
    uniquecondi=unique(conditiontmp);
    lengthcondi=length(uniquecondi);
    clength=ceil(sqrt(lengthcondi));
    rlength=floor(sqrt(lengthcondi));
    while clength*rlength<lengthcondi
        clength=clength+1;
    end
    if lengthcondi>maxsub
        k=1;
        clength=5;
        rlength=4;
        for isub=1:lengthcondi
            figure(h(k))
            i=mod(isub,maxsub);
            if mod(isub,maxsub)==0
                i=maxsub;
                k=k+1;
                h(k)=figure('NumberTitle','off','Name',['Mean fixation bias - ' CName{ic} ' continue']);
                figure(h(k-1))
            end
            subplot(rlength,clength,i)
            idxtmp=strcmp(conditiontmp,uniquecondi(isub));
            TempMap=squeeze(nanmean(FixMap(idxtmp,:,:),1));
            imagesc(TempMap)
            title(uniquecondi{isub})
            axis off,axis equal
        end
    else
        for isub=1:lengthcondi
            subplot(rlength,clength,isub)
            idxtmp=strcmp(conditiontmp,uniquecondi(isub));
            TempMap=squeeze(nanmean(FixMap(idxtmp,:,:),1));
            imagesc(TempMap)
            title(uniquecondi{isub})
            axis off,axis equal
        end
    end
end
%% text output
label=cell(length(MeasureM),length(CName)-1);
strlabel=cell(length(MeasureM),1);
for icc=1:length(CName)-1
    label(:,icc)=cellstr(eval(['ConditionM.' CName{icc}]));
end
for ii=1:length(label)
    strlabel{ii,:}=strjoin(label(ii,:),'_');
end
% [CatePredictor,btmp]=unique(strlabel,'rows');

%%
%%
h(1)=figure('NumberTitle','off','Name',['Mean fixation bias - JointCondition']);
conditiontmp=strlabel;
uniquecondi=unique(conditiontmp);
lengthcondi=length(uniquecondi);
clength=ceil(sqrt(lengthcondi));
rlength=floor(sqrt(lengthcondi));
while clength*rlength<lengthcondi
    clength=clength+1;
end
if lengthcondi>maxsub
    k=1;
    clength=5;
    rlength=4;
    for isub=1:lengthcondi
        figure(h(k))
        i=mod(isub,maxsub);
        if mod(isub,maxsub)==0
            i=maxsub;
            k=k+1;
            h(k)=figure('NumberTitle','off','Name',['Mean fixation bias - JointCondion continue']);
            figure(h(k-1))
        end
        subplot(rlength,clength,i)
        idxtmp=strcmp(conditiontmp,uniquecondi(isub));
        TempMap=squeeze(nanmean(FixMap(idxtmp,:,:),1));
        imagesc(TempMap)
        title(uniquecondi{isub})
        axis off,axis equal
    end
else
    for isub=1:lengthcondi
        subplot(rlength,clength,isub)
        idxtmp=strcmp(conditiontmp,uniquecondi(isub));
        TempMap=squeeze(nanmean(FixMap(idxtmp,:,:),1));
        imagesc(TempMap)
        title(uniquecondi{isub})
        axis off,axis equal
    end
end

