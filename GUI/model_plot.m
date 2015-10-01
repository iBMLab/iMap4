function model_plot(LMMmap,FixMap,pathname)
global StatMap StatMap_c
screensize = get(0,'ScreenSize');
size_figure = [400,100];
xpos = ceil((screensize(3)-size_figure(2))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-size_figure(1))/2); % center the figure on the screen vertically

Mask=squeeze(~isnan(LMMmap.MSE));
S.fh = figure('units','pixels',...
    'position',[xpos ypos size_figure(1) size_figure(2)],...
    'menubar','none',...
    'name','Model Statistics',...
    'numbertitle','off',...
    'resize','off','ToolBar','none');
t = uitoolbar(S.fh);
% Read an image
[img,map] = imread(fullfile(matlabroot,'toolbox','matlab','icons','helpicon.gif'));
icon = ind2rgb(img,map);

uitoggletool(t,'CData',icon,'TooltipString','Help','ClickedCallback',@help_button);


set(S.fh,'color','w'); % white background

sz_pshb1 = [85,30];

%pshb_model
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Model Fitting','position',[20 50 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',@model);

%pshb_anova
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','ANOVA','position',[115 50 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',@anova);

%pshb_contrast
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Linear Contrast','position',[210 50 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',@contrast);

%pshb_posthoc
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Post-Hoc','position',[305 50 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',@post_hoc);

%pshb_multiple_comparison_contrast
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Multiple Comparison Correction','position',[30 10 160 sz_pshb1(2)],...
    'Callback',@multiple_comparison_correction);


%pshb_display_results
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Display Results','position',[220 10 160 sz_pshb1(2)],...
    'Callback',@display_results);

uiwait(gcf)

    function model(~,~)
        
        %plot model fitting
        opt1.type='model';
        
        %add the if StatMap.mccopt and the StatMap.opt.type
        
        %perform contrast
        StatMap=imapLMMcontrast(LMMmap,opt1);
        save(strcat(pathname,'StatMap_model'),'StatMap','-v7.3')
        
        %option for LMMdisplay
        %normalized,backgroundfile,colourmap,colormaprange,distplot
        prompt={'background file','colormap','map range for optimal display (0 to 1)','map Value distribution (1 to display)'};
        name=' Options for Display';
        numlines = [1 50;1 50;1 50;1 50];
        defaultanswer={'','''jet''','',''};
        options=inputdlg(prompt,name,numlines,defaultanswer);
        if isempty(options)==0
            normalized      = 0;
            backgroundfile  = char(options{1});
            cmap            = options{2};
            colormaprange   = str2double(options{3});
            distplot        = str2double(options{4});
        else
            normalized      = 0;
            backgroundfile  = char(defaultanswer{1});
            cmap            = defaultanswer{2};
            colormaprange   = str2double(defaultanswer{3});
            distplot        = str2double(defaultanswer{4});
        end
        %output figure;
        % imapLMMdisplay(StatMap,0)
        imapLMMdisplay(StatMap,normalized,backgroundfile(2:end-1),cmap(2:end-1),colormaprange,distplot);
    end

    function anova(~,~)
        %option of imapLMMcontrast
        prompt={'alpha threshold'};
        name='Please select the p value';
        numlines = [1 50];
        defaultanswer={'0.05'};
        options=inputdlg(prompt,name,numlines,defaultanswer);
        if isempty(options)==0
            opt.type     = 'fixed';
            opt.alpha    = str2double(options{1});
        else
            opt.type     = 'fixed';
            opt.alpha    = str2double(defaultanswer{1});
        end
        % perform contrast
        StatMap=imapLMMcontrast(LMMmap,opt);
        
        %ouput warning
        h = warndlg(' The displayed maps are the ANOVA of the fixed effect before multiple comparison correction, the correction on the statistic map will follow shortly');
        set(h,'color','white');
        uiwait(gcf)
        % save StatMap as StatMap_ANOVA.mat
        save(strcat(pathname,'StatMap_ANOVA'),'StatMap','-v7.3')
        
        %         %option for LMMdisplay
        %         %normalized,backgroundfile,colourmap,colormaprange,distplot
        %         prompt={'Normalize map values (1 to normalize all maps)','Background file','Colormap (as being used in Matlab)','map range for optimal display (0 to 1)','map Value distribution (1 to display)'};
        %         name='Options for Display';
        %         numlines = [1 50;1 50;1 50;1 50;1 50];
        %         defaultanswer={'','','''jet''','',''};
        %         options=inputdlg(prompt,name,numlines,defaultanswer);
        %         normalized      = str2double(options{1});
        %         cmap            = options{3};
        %         backgroundfile  = char(options{2});
        %         colormaprange   = str2double(options{4});
        %         distplot        = str2double(options{5});
        
        %output figure;
        imapLMMdisplay(StatMap,0)
        % imapLMMdisplay(StatMap,normalized,backgroundfile(2:end-1),cmap(2:end-1),colormaprange,distplot);
        
        % go to the Multiple comparison tab
        StatMap_c=multiple_comparison_tab(LMMmap,FixMap,StatMap,pathname);
        
    end

    function contrast(~,~,~)
        if isfield(LMMmap,'SinglePred')
            size_figure = [300,300];
            screensize = get(0,'ScreenSize');
            sz = [185,180];
            %sz = [300,200];
            xpos = ceil((screensize(3)-sz(2))/2); % center the figure on the screen horizontally
            ypos = ceil((screensize(4)-sz(1))/2); % center the figure on the screen vertically
            
            S.fh = figure('units','pixels','position',[xpos, ypos size_figure(1) size_figure(2)],...
                'menubar','none',...
                'name','Contrast Selection',...
                'numbertitle','off',...
                'resize','off');
            %movegui(S.fh,'center')
            set(gcf,'color','w');
            
            r = {'Kindly select the predictors you would like to perform linear contrast'};
            S.text_general = uicontrol('parent',S.fh,'style','text','string',r,'position',[50 310  500 30],...
                'HorizontalAlignment','center','FontWeight','Bold','Backgroundcolor','white');
            
            n = length(LMMmap.SinglePred.CatePredictor);
            levels = [LMMmap.SinglePred.CatePredictor, num2cell(false(n,1)),num2cell(false(n,1))];
            %sz_table_random = [160,150];
            columnname = {'Categorical|Predictors','Positive|contrast','Negative|contrast'};
            column_type = cell(n,1);
            column_type(:) = cellstr('logical');
            columnformat = ['char',column_type'];
            
            
            categorical_predictor_table = uitable(S.fh,'position',[30,30,226,250],'Units','normalized',...
                'data', levels,...
                'ColumnEditable', [false(1,1),true(n,1)'],...
                'columnFormat', columnformat,...
                'rowName',[],'RowStriping','off',...
                'ColumnName',columnname);
            
            sz_pshb1 = [70,20];
            
            uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
                'String','Validate','position',[60 8 sz_pshb1(1) sz_pshb1(2)],...
                'Callback',{@validate_contrast});
            
            % create push button 'validate'
            uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
                'String','Cancel','position',[150 8 70 20],...
                'Callback',@cancel_contrast);
            uiwait(gcf)
        else
            errordlg('You need to estimate the predictor beta in the imapLMM')
            uiwait(gcf)
            return
        end
        
        function validate_contrast(~,~)
            uiresume(gcf)
            opt=struct;% clear structure
            
            retreive_data = get(categorical_predictor_table,'data');
            if ishandle(S)
                delete(S.fh)
            end
            positive_contrast = retreive_data(:,2);
            negative_contrast = retreive_data(:,3);
            
            n_positive = find(cell2mat(positive_contrast)==1);
            n_negative = find(cell2mat(negative_contrast)==1);
            
            if isempty(n_negative)==0 && isempty(n_positive)==0 % positive and negative contrast
                
                opt.c = zeros(1,size(retreive_data,1));
                opt.c(n_positive) = 1;
                opt.c(n_negative) = -1;
                %name = strjoin(retreive_data([n_positive;n_negative],1)','/');
                %opt.name = name ;
                if sum(opt.c)~=0 % balance the vector of contrast
                    newc=opt.c;
                    %errordlg('non-balanced vector!');
                    check = sign(sum(newc));
                    if check<0
                        k = find(newc>0);
                    else
                        k = find(newc<0);
                    end
                    opt.c(k) = (-check).*(abs(newc(k)*abs(sum(newc))/length(k)) + abs(newc(k)));
                end
                
                %option of imapLMMcontrast
                prompt={'h','alpha','Name of the contrast'};
                dlg_title = 'Options for imapLMMcontrast';
                num_lines = [1 60];
                defaultanswer={'0','0.05','Contrast'};
                options = inputdlg(prompt,dlg_title,num_lines,defaultanswer);
                opt.type     = 'predictor beta';
                opt.h       =  str2double(options{1});
                opt.alpha    = str2double(options{2});
                if ~isempty(options{3})
                    opt.name = options(3);
                end
                if isempty(options{1})
                    opt.h=0;
                end
            else
                if length(n_negative)==1 || length(n_positive)==1 % one negative/positive contrast
                    
                    opt.c = zeros(1,size(retreive_data,1));
                    opt.c(n_positive) = 1;
                    opt.c(n_negative) = 1;
                    opt.h = mean(mean(FixMap(:,Mask)));
                    if length(n_positive)==1
                        name = retreive_data(cell2mat(positive_contrast)==1,1);
                        opt.name = name;
                        opt =  function_contrast_one_tail(opt);
                        if isempty(opt.name)
                            opt.name = name;
                        end
                    else
                        name = retreive_data(cell2mat(negative_contrast)==1,1);
                        opt.name = name;
                        opt =  function_contrast_one_tail(opt);
                        if isempty(opt.name)
                            opt.name = name;
                        end
                    end
                    opt.onetail='>';
                    %msgbox(strcat('We will perform a one tail t.test of predictors',num2str(''))
                    
                else
                    if isempty(n_negative)==1 && isempty(n_positive)==1
                        h = errordlg('You did not select any predictor');
                        set(h,'color','white')
                        uiwait(gcf)
                    else % multiple positive or negative contrast
                        h2 = warndlg('Non-balanced vector');
                        set(h2,'color','white')
                        uiwait(gcf)
                        opt.c = zeros(1,size(retreive_data,1));
                        opt.c(n_positive) = 1;
                        opt.c(n_negative) = -1;
                        
                        %option of imapLMMcontrast
                        prompt={'h','alpha','Name of the contrast'};
                        dlg_title = 'Options for imapLMMcontrast';
                        num_lines = [1 60];
                        defaultanswer={'0','0.05',''};
                        options = inputdlg(prompt,dlg_title,num_lines,defaultanswer);
                        opt.type     = 'predictor beta';
                        if isempty(options)==0
                            opt.h       =  str2double(options{1});
                            opt.alpha    = str2double(options{2});
                            if ~isempty(options{3})
                                opt.name = options(3);
                            end
                            if isempty(options{1})
                                opt.h=0;
                            end
                        else
                            opt.h        =  str2double(defaultanswer{1});
                            opt.alpha    =  str2double(defaultanswer{2});
                        end
                        
                    end
                end
            end
            if ishandle(S.fh)
                delete(S.fh)
            end
            if isfield(opt,'type')==1
                % perform contrast
                [StatMap]=imapLMMcontrast(LMMmap,opt);
                
                %ouput warning
                h = warndlg('The displayed maps are the Linear Contrast of the fixed effect, please process to perform multiple comparison correction');
                set(h,'color','white');
                uiwait(gcf)
                %save StatMap as StatMap_.mat asked user to rename
                prompt = {'Please rename the StatMap',};
                dlg_title = 'Saving StatMap';
                num_lines = 1;
                def = {'StatMap_Contrast_'};
                answer = [];
                while isempty(answer)
                    answer = inputdlg(prompt,dlg_title,num_lines,def);
                end
                
                save(strcat(pathname,answer{1}), 'StatMap','-v7.3');
                
                %output figure;
                imapLMMdisplay(StatMap,0)
                % imapLMMdisplay(StatMap,normalized,backgroundfile(2:end-1),cmap(2:end-1),colormaprange,distplot);
                
                StatMap_c=multiple_comparison_tab(LMMmap,FixMap,StatMap,pathname);
            end
        end
        
        function cancel_contrast(~,~)
            delete(gcf)
        end
    end

    function post_hoc(~,~)
        if isempty(StatMap_c)==0
            
            %load RawMap
            [filename, pathname] = uigetfile('*.mat',' Select RawMap to proceed');
            fullpathname = strcat(pathname,filename); %user selects only one file
            if isempty(fullpathname)
            h = errordlg('No file selected');
            set(h,'color','white');
            else
            RawmapMAT=load(strcat(pathname,filename));
            try
                [Posthoc]=imapLMMposthoc(StatMap_c,RawmapMAT.RawMap,LMMmap,'mean',1);
            catch
                errordlg('Please select the right file')
            end
            
            %save StatMap as StatMap_.mat asked user to rename
            prompt = {'Please rename the Posthoc',};
            dlg_title = 'Saving Posthoc';
            num_lines = 1;
            def = {'Posthoc_'};
            answer = [];
            while isempty(answer)
                answer = inputdlg(prompt,dlg_title,num_lines,def);
            end
            
            save(strcat(pathname,answer{1}), 'Posthoc','-v7.3');
            end
        else if isempty(StatMap)==0
                
                %load RawMap
                [filename, pathname] = uigetfile('*.mat',' Select RawMap to proceed');
                RawmapMAT=load(strcat(pathname,filename));
                try
                    [Posthoc]=imapLMMposthoc(StatMap_c,RawmapMAT.RawMap,LMMmap,'mean',1);
                catch
                    errordlg('Please select the right file')
                end
                
                %save StatMap as StatMap_.mat asked user to rename
                prompt = {'Please rename the Posthoc',};
                dlg_title = 'Saving Posthoc';
                num_lines = 1;
                def = {'Posthoc_'};
                answer = [];
                while isempty(answer)
                    answer = inputdlg(prompt,dlg_title,num_lines,def);
                end
                
                save(strcat(pathname,answer{1}), 'Posthoc','-v7.3');
            else
              h =   errordlg('You need to perform ANOVA or Linear Contrast before doing the Post-Hoc');
             set(h,'color','white');
            end
        end
    end

    function display_results(~,~)
        
        [filename, pathname] = uigetfile('*.mat',' Select StatMap to proceed');
        fullpathname = strcat(pathname,filename); %user selects only one file
        
        if isempty(fullpathname)
            h = errordlg('No file selected');
            set(h,'color','white');
        else
            loadMAP=load(strcat(pathname,filename));
            if isfield(loadMAP,'StatMap')
                StatMaptoplot=loadMAP.StatMap;
            elseif isfield(loadMAP,'StatMap_c')
                StatMaptoplot=loadMAP.StatMap_c;
            else
                h = errordlg('The StatMap is mandatory to proceed');
                set(h,'color','white');
                return
            end
            if strcmp(StatMaptoplot.opt.type,'model')
                %option for LMMdisplay
                %normalized,backgroundfile,colourmap,colormaprange,distplot
                prompt={'background file','colormap','map range for optimal display (0 to 1)','map Value distribution (1 to display)'};
                name=' Options for Display';
                numlines = [1 50;1 50;1 50;1 50];
                defaultanswer={'','''jet''','',''};
                options=inputdlg(prompt,name,numlines,defaultanswer);
                normalized      = 0;
                backgroundfile  = char(options{1});
                cmap            = options{2};
                colormaprange   = str2double(options{3});
                distplot        = str2double(options{4});
            else
                %option for LMMdisplay
                %normalized,backgroundfile,colourmap,colormaprange,distplot
                prompt={'Normalize map values (1 to normalize all maps)','Background file','Colormap (as being used in Matlab)','map range for optimal display (0 to 1)','map Value distribution (1 to display)'};
                name='Options for Display';
                numlines = [1 50;1 50;1 50;1 50;1 50];
                defaultanswer={'','','''jet''','',''};
                options=inputdlg(prompt,name,numlines,defaultanswer);
                if isempty(options)==0
                normalized      = str2double(options{1});
                cmap            = options{3};
                backgroundfile  = char(options{2});
                colormaprange   = str2double(options{4});
                distplot        = str2double(options{5});
                else
                normalized      = str2double(defaultanswer{1});
                cmap            = defaultanswer{3};
                backgroundfile  = char(defaultanswer{2});
                colormaprange   = str2double(defaultanswer{4});
                distplot        = str2double(defaultanswer{5});    
                end
            end
            %output figure;
            % imapLMMdisplay(StatMap,0)
            imapLMMdisplay(StatMaptoplot,normalized,backgroundfile(2:end-1),cmap(2:end-1),colormaprange,distplot);
        end
        
        
    end
    function multiple_comparison_correction(~,~)
        
        [filename, pathname] = uigetfile('*.mat',' Select StatMap to proceed');
        fullpathname = strcat(pathname,filename);
        if isempty(fullpathname)
            errordlg('No file selected')
        else
            load_StatMap = load(strcat(pathname,filename));
            if isfield(load_StatMap,'StatMap')==0
                h = errordlg('The StatMap is mandatory to proceed');
                set(h,'color','white')
                uiwait(gcf)
            else
                StatMap_c=multiple_comparison_tab(LMMmap,FixMap,load_StatMap.StatMap,pathname);
            end
        end
        
        
    end
end

function help_button(~,~)
if ispc
    winopen('manual_GUI.pdf');
elseif ismac
    system(['open manual_GUI.pdf']);
else
    errordlg('Please open the Guidebook manually in ./matlab/Apps/iMAP/GUI')
end
end