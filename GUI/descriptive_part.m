function descriptive_part(varargin)
scriptName = mfilename('fullpath');
[currentpath, ~, ~]= fileparts(scriptName);
addpath(currentpath);
%
if nargin==0
    [filename, pathname] = uigetfile('*.mat',' Select DescriptvM, FixMap to proceed','MultiSelect','on');
    if iscell(filename)==0 %user selects only multiple file
        h = errordlg('DescriptvM and FixMap files need to be selected');
        %set(h,'color','white');
        uiwait(h)
    else
        % check of the selected files
        h = waitbar(0,'Please wait loading files is in progress','Name','Loading files','color','w');
        
        for i = 1:length(filename)
            load(strcat(pathname,filename{i}))
            waitbar(i/length(filename))
        end
        if ishandle(h)
            delete(h)
        end
    end
elseif nargin==2
    DescriptvM = varargin{1};
    FixMap     = varargin{2};
elseif nargin==3
    DescriptvM = varargin{1};
    FixMap     = varargin{2};
    opt        = varargin{3};
end
%%
if exist('DescriptvM','var') && exist('FixMap','var')
    % close all
    % indx=find(DescriptvM.FixNum>200);
    % FixMap(indx,:,:)=[];
    % RawMap(indx,:,:)=[];
    
    % DescriptvM(indx,:)=[];
    % PredictorM(indx,:)=[];
    %% Create a folder to save all figures
    if exist('Descriptive_STAT', 'file')==0
        mkdir('Descriptive_STAT');
        addpath('Descriptive_STAT');
    end
    %%
    ConditionM=DescriptvM(:,6:end);
    MeasureM=DescriptvM(:,1:5);
    nbins=50;
    % overall distribution
    
    screensize = get(0,'ScreenSize');
    xpos = ceil((screensize(3)-screensize(3))/2); % center the figure on the screen horizontally
    ypos = ceil((screensize(4)-screensize(4))/2); % center the figure on the screen vertically
    
    h1 = figure('NumberTitle','off','Name','Eye movement measurement distribution',...
        'position',[xpos,ypos,screensize(3)/2,screensize(4)]);
    % given that fixation number could be consider as a possion process, the
    % parameters could be fit with a Gamma distribution. Fixation
    % duration/path length could both be consider some how as the waiting time
    % of the possion process, thus should be appropriate to fit also a gamma
    % distribution
    subplot(3,2,1)
    histfit(MeasureM.FixNum,nbins,'gamma');
    title('Fixation Number')
    
    subplot(3,2,3)
    histfit(MeasureM.sumFixDur,nbins,'gamma');
    title('Fixation Duration(SUM)')
    
    subplot(3,2,4)
    histfit(MeasureM.meanFixDur,nbins,'gamma');
    title('Fixation Duration(MEAN)')
    
    subplot(3,2,5)
    histfit(MeasureM.totalPathLength,nbins,'gamma');
    title('Path Length(TOTAL)')
    
    subplot(3,2,6)
    histfit(MeasureM.meanPathLength,nbins,'gamma');
    title('Path Length(MEAN)')
    if isunix
        print(strcat('Descriptive_STAT/Eye_Mvt_Measure_Dist'),'-dpng')
    else
        print(strcat('Descriptive_STAT\Eye_Mvt_Measure_Dist'),'-dpng')
    end
    
    % Mean of Fixation Map
    h2 = figure('NumberTitle','off','Name','Mean Fixation Intensity',...
        'position',[xpos,ypos,screensize(3),screensize(4)]);
    title('Mean Fixation')
    FixationMean=squeeze(nanmean(FixMap(:,:,:)));
    imagesc(FixationMean);
    axis off,axis equal
    if isunix
        print(strcat('Descriptive_STAT/Mean_Intensity_Map'),'-dpng')
    else
        print(strcat('Descriptive_STAT\Mean_Intensity_Map'),'-dpng')
    end
    close(h1,h2)
    
    % text output
    txtimapout(MeasureM,[],[],1)
    %% box plot for each catigorical condition
    
    CName = ConditionM.Properties.VarNames;
    
    for ic=1:length(CName)
        if isunix
            mkdir(strcat('Descriptive_STAT/',CName{ic}));
            addpath(strcat('Descriptive_STAT/',CName{ic}));
        else
            mkdir(strcat('Descriptive_STAT\',CName{ic}));
            addpath(strcat('Descriptive_STAT\',CName{ic}));
        end
        h1 = figure('NumberTitle','off','Name',['Condition ' CName{ic}],...
            'position',[xpos,ypos,screensize(3)/2,screensize(4)]);
        subplot(3,2,1)
        Conditiontmp=eval(['ConditionM.' CName{ic}]);
        boxplot(MeasureM.FixNum,Conditiontmp)
        title('Fixation Number')
        
        subplot(3,2,3)
        boxplot(MeasureM.sumFixDur,Conditiontmp)
        title('Fixation Duration(SUM)')
        
        subplot(3,2,4)
        boxplot(MeasureM.meanFixDur,Conditiontmp)
        title('Fixation Duration(MEAN)')
        
        subplot(3,2,5)
        boxplot(MeasureM.totalPathLength,Conditiontmp)
        title('Path Length(TOTAL)')
        
        subplot(3,2,6)
        boxplot(MeasureM.meanPathLength,Conditiontmp)
        title('Path Length(MEAN)')
        if isunix
            print(strcat('Descriptive_STAT/',CName{ic},'/',CName{ic},'_Boxplot'),'-dpng')
        else
            print(strcat('Descriptive_STAT\',CName{ic},'\',CName{ic},'_Boxplot'),'-dpng')
        end
        close(h1)
        
        txtimapout(MeasureM,Conditiontmp,CName{ic},2)
    end
    
    %% mean fixation map
    maxsub=20;
    for ic=1:length(CName)-1
        %%
        h(1)=figure('NumberTitle','off','Name',['Mean fixation bias - ' CName{ic}],...
            'position',[xpos,ypos,screensize(3),screensize(4)]);
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
            pp=0;
            clength=5;
            rlength=4;
            for isub=1:lengthcondi
                figure(h(k))
                i=mod(isub,maxsub);
                if mod(isub,maxsub)==0;i=maxsub;k=k+1;pp=1;end
                subplot(rlength,clength,i)
                try
                    idxtmp=conditiontmp==uniquecondi(isub);
                catch
                    idxtmp=strcmp(conditiontmp,uniquecondi(isub));
                end
                TempMap=squeeze(nanmean(FixMap(idxtmp,:,:),1));
                imagesc(TempMap)
                title(char(uniquecondi(isub)))
                axis off,axis equal
                
                if pp==1
                    if isunix
                        print(strcat('Descriptive_STAT/',CName{ic},'/',CName{ic},'_Mean_Intensity_Map',num2str(k-1)),'-dpng')
                    else
                        print(strcat('Descriptive_STAT\',CName{ic},'\',CName{ic},'_Mean_Intensity_Map',num2str(k-1)),'-dpng')
                    end
                    h(k)=figure('NumberTitle','off','Name',['Mean_FixaionIntensity_Map - ' CName{ic} ' continue'],...
                        'position',[xpos,ypos,screensize(3),screensize(4)]);
                    close(h(k-1))
                    pp=0;
                end
            end
            if isunix
                print(strcat('Descriptive_STAT/',CName{ic},'/',CName{ic},'_Mean_Intensity_Map',num2str(k)),'-dpng')
            else
                print(strcat('Descriptive_STAT\',CName{ic},'\',CName{ic},'_Mean_Intensity_Map',num2str(k)),'-dpng')
            end
            close(gcf)
        else
            for isub=1:lengthcondi
                subplot(rlength,clength,isub)
                try
                    idxtmp=conditiontmp==uniquecondi(isub);
                catch
                    idxtmp=strcmp(conditiontmp,uniquecondi(isub));
                end
                TempMap=squeeze(nanmean(FixMap(idxtmp,:,:),1));
                imagesc(TempMap)
                title(char(uniquecondi(isub)))
                axis off,axis equal
            end
            % save figures
            if isunix
                print(strcat('Descriptive_STAT/',CName{ic},'/',CName{ic},'_Mean_Intensity_Map'),'-dpng')
            else
                print(strcat('Descriptive_STAT\',CName{ic},'\',CName{ic},'_Mean_Intensity_Map'),'-dpng')
            end
            %close figure
            close(gcf);
        end
    end
    CName2 = 0; % initialize CName2
    while isempty(CName2)==0
        selection_conditions = 1;
        while (selection_conditions)==1
            [CName2,selection_conditions] = function_select_predictors_stat(ConditionM);
        end
        if isempty(CName2)==0
            CNameall=strjoin(CName2','-');
            if isunix
                mkdir('Descriptive_STAT/',CNameall);
                addpath('Descriptive_STAT/',CNameall);
            else
                mkdir('Descriptive_STAT\',CNameall);
                addpath('Descriptive_STAT\',CNameall);
            end
            %%
            label=cell(length(MeasureM),length(CName2));
            strlabel=cell(length(MeasureM),1);
            for icc=1:length(CName2)
                try
                    label(:,icc)=cellstr(eval(['ConditionM.' CName2{icc}]));
                catch
                    for ii=1:length(MeasureM)
                        label{ii,icc}=num2str(eval(['ConditionM.' CName2{icc} '(ii)']));
                    end
                end
            end
            for ii=1:length(label)
                strlabel{ii,:}=strjoin(label(ii,:),'-');
            end
            % [CatePredictor,btmp]=unique(strlabel,'rows');
            
            %% Joint Condition
            h(1)=figure('NumberTitle','off','Name',['Mean Fixation Intensity - JointCondition'],...
                'position',[xpos,ypos,screensize(3),screensize(4)]);
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
                pp=0;
                clength=5;
                rlength=4;
                for isub=1:lengthcondi
                    figure(h(k))
                    i=mod(isub,maxsub);
                    if mod(isub,maxsub)==0;i=maxsub;k=k+1;pp=1;end
                    subplot(rlength,clength,i)
                    idxtmp=strcmp(conditiontmp,uniquecondi(isub));
                    TempMap=squeeze(nanmean(FixMap(idxtmp,:,:),1));
                    imagesc(TempMap)
                    title(uniquecondi{isub})
                    axis off,axis equal
                    
                    if pp==1
                        if isunix
                            print(strcat('Descriptive_STAT/',CNameall,'/',CNameall,'_Mean_Intensity_Map',num2str(k-1)),'-dpng')
                        else
                            print(strcat('Descriptive_STAT\',CNameall,'\',CNameall,'_Mean_Intensity_Map',num2str(k-1)),'-dpng')
                        end
                        h(k)=figure('NumberTitle','off','Name',['Mean Fixation Intensity - JointCondion continue'],...
                            'position',[xpos,ypos,screensize(3),screensize(4)]);
                        close(h(k-1))
                        pp=0;
                    end
                end
                if isunix
                    print(strcat('Descriptive_STAT/',CNameall,'/',CNameall,'_Mean_Intensity_Map',num2str(k)),'-dpng')
                else
                    print(strcat('Descriptive_STAT\',CNameall,'\',CNameall,'_Mean_Intensity_Map',num2str(k)),'-dpng')
                end
                close(gcf)
            else
                for isub=1:lengthcondi
                    subplot(rlength,clength,isub)
                    idxtmp=strcmp(conditiontmp,uniquecondi(isub));
                    TempMap=squeeze(nanmean(FixMap(idxtmp,:,:),1));
                    imagesc(TempMap)
                    title(uniquecondi{isub})
                    axis off,axis equal
                end
                if isunix
                    print(strcat('Descriptive_STAT/',CNameall,'/',CNameall,'_Mean_Intensity_Map'),'-dpng')
                else
                    print(strcat('Descriptive_STAT\',CNameall,'\',CNameall,'_Mean_Intensity_Map'),'-dpng')
                end
                close(gcf);
            end
            % %% savefigures
            % hfigs = get(0, 'children') ;
            % hfigs = sort(hfigs);%Get list of figures
            % %hfigs(hfigs ==handles.figure1)=[];
            % close(hfigs(1:end-1)) % keep the figure of the linear model open
            
            %% box plot for  joint condition
            h(1) = figure('NumberTitle','off','Name','Eye movement measurement Boxplot',...
                'position',[xpos,ypos,screensize(3)/2,screensize(4)]);
            
            subplot(3,2,1)
            boxplot(MeasureM.FixNum,conditiontmp);
            title('Fixation Number')
            
            subplot(3,2,3)
            boxplot(MeasureM.sumFixDur,conditiontmp);
            title('Fixation Duration(SUM)')
            
            subplot(3,2,4)
            boxplot(MeasureM.meanFixDur,conditiontmp);
            title('Fixation Duration(MEAN)')
            
            subplot(3,2,5)
            boxplot(MeasureM.totalPathLength,conditiontmp);
            title('Path Length(TOTAL)')
            
            subplot(3,2,6)
            boxplot(MeasureM.meanPathLength,conditiontmp);
            title('Path Length(MEAN)')
            if isunix
                print(strcat('Descriptive_STAT/',CNameall,'/',CNameall,'_Boxplot'),'-dpng');
            else
                print(strcat('Descriptive_STAT\',CNameall,'\',CNameall,'_Boxplot'),'-dpng');
            end
            close(gcf)
            
            h1= msgbox('All Figures are closed and saved in the folder Descriptive_STAT');
            uiwait(h1)
            
            txtimapout(MeasureM,conditiontmp,CNameall,2)
        end
    end
    %% Modelling Analysisclear all
    if ~exist('opt')
        modelling_part()
    end
else
    errordlg('Please select the correct Files: DescriptivM and FixMap')
end