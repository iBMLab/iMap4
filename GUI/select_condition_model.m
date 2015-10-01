function [formula,opt,opttext] = select_condition_model(PredictorM,file)
formula = []; opt.singlepredi=1;opttext={'DummyVarCoding','effect','FitMethod','REML'};
idx_formula = 0; 
size_figure = [600,350];
screensize = get(0,'ScreenSize');
sz = [185,180];
%sz = [300,200];
xpos = ceil((screensize(3)-sz(2))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-sz(1))/2); % center the figure on the screen vertically

S.fh = figure('units','pixels','position',[xpos, ypos size_figure(1) size_figure(2)],...
    'menubar','none',...
    'name','Selection of fixed and random Predictors',...
    'numbertitle','off',...
    'resize','off');
%movegui(S.fh,'center')
set(gcf,'color','w');

r = {'Kindly select the predictors you would like to include as fixed/random effect in the LMM'};
S.text_general = uicontrol('parent',S.fh,'style','text','string',r,'position',[20 310  500 30],...
    'HorizontalAlignment','center','FontWeight','Bold','Backgroundcolor','white');

%random predictors
random_predictors = [PredictorM.Properties.VarNames];
random_predictors = fliplr(random_predictors);
% duration = Pr(end-1);%get the duration
% sbjST    = condition_name(end);%get the subject

%%%random predictors
n = length(random_predictors);
levels = [random_predictors', [num2cell(true(1,1));num2cell(false(n-1,1))]];
%sz_table_random = [160,150];
columnname = {' Random Predictors ','CheckBox'};
column_type = cell(n,1);
column_type(:) = cellstr('logical');
columnformat = ['char',column_type'];


random_predictors_table = uitable(S.fh,'position',[310,60,190,250],'Units','normalized',...
    'data', levels,...
    'ColumnEditable', [false(1,1),true(n,1)'],...
    'columnFormat', columnformat,...
    'rowName',[],'RowStriping','off',...
    'ColumnName',columnname,'ColumnWidth', 'auto');

%create interaction conditions
interaction_condition = [];
if isempty(strfind(file,'single_trial'))==0 %single trial FixMap
random_predictors(1:2) =  [];% first two variables are subject and duration
random_predictors(end) =  [];% last variable is the Trial (since we fliped!)
else
random_predictors(1:2) =  [];% first two variables are subject and duration
    
end
random_predictors = fliplr(random_predictors); %kind of sorting

for i=1:length(random_predictors)
    if i==1
        combination1 = combnk(random_predictors,i);
        if size(combination1,2)> size(combination1,1)
            combination1 = combination1';
        end
    else
        comb = combnk(random_predictors,i);
        comb = flipud(comb);
        combination = cell(size(comb,1),1);
        for jj=1:size(comb,1)
            combination(jj,1) = cellstr(implode(comb(jj,:),':'));
        end
        interaction_condition = [interaction_condition;combination];
        
    end
end
interaction_condition = [combination1;interaction_condition];

condition_name =   interaction_condition;
n_condition_name = length(condition_name);
levels = [condition_name,num2cell(true(n_condition_name,1))];
columnname = {'Fixed Predictors','CheckBox'};
column_type = cell(n_condition_name,1);
column_type(:) = cellstr('logical');
columnformat = ['char',column_type'];

fixed_predictors_table = uitable('position',[30,60,188,250],'Units','normalized',...
    'data', levels,...
    'ColumnEditable', [false(1,1),true(n_condition_name,1)'],...
    'columnFormat', columnformat,...
    'rowName',[],'RowStriping','off',...
    'ColumnName',columnname,'ColumnWidth','auto');

% set(handles.figure1, 'pointer', 'arrow')

% create push button 'uncheck all'
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Uncheck all','position',[40 15 80 20],...
    'Callback',{@call_setall,0});

% create push button 'check all'
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Check all','position',[130 15 80 20],...
    'Callback',{@call_setall,1});

uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Edit Formula','position',[270 15 80 20],...
    'Callback',@see_formula);

uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Run Model','position',[450 15 80 20],...
    'Callback',{@run_model});

%  %create push button 'cancel'
%  uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
%          'String','Cancel','position',[450 15 80 20],...
%          'Callback',@cancel_model);

% create push button 'option'
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Option Model','position',[360 15 80 20],...
    'Callback',@option_model);

uiwait(gcf)

    function run_model(~,~)
        uiresume(gcf)
        if idx_formula == 0
        retrieve_random_table = get(random_predictors_table,'data');
        retrieve_fixed_table = get(fixed_predictors_table,'data');
        random_factors = retrieve_random_table(cell2mat(retrieve_random_table(:,2))==1,1);
        fixed_factors = retrieve_fixed_table(cell2mat(retrieve_fixed_table(:,2))==1,1);
        
        factors = [random_factors;fixed_factors];
        % % change table to nominal for modelling
        
        for h = 1:length(factors)
            
            if isempty(strfind(factors{h},':'))==1
                idx = unique(PredictorM.(factors{h}));
                if length(idx) ==1
                    h1 = warndlg(strcat('The factor ',{' '}, factors{h},{' '},' has been excluded because it has only one level'));
                    set(h1,'color','white')
                    fixed_factors(strncmp(factors{h}, fixed_factors, length(factors{h})))=[];
                    %else
                    %PredictorM.(factors{h}) = nominal(PredictorM.(factors{h}));
                end
            end
        end
        
        random_part = cell(1,length(random_factors));
        for j1=1:length(random_factors)
%             if j1==length(random_factors)
                random_part{1,j1} = strcat(random_factors{j1},')');
%             else
%                 random_part{1,j1} = strcat(random_factors{j1},')+(1|');
%             end
        end
        random_part = strcat('+(1|',random_part);
        fixed_part = implode(fixed_factors,'+');
        formula = strcat('PixelIntensity ~ ',' ',fixed_part,strjoin(random_part));
        end
        if ishandle(S.fh)
            delete(S.fh)
        end
    end

    function call_setall(~,~,value)
        % check (value=1) or uncheck (value=0) all checkboxes
        if value ==1
            fixed_predictors_table = uitable('position',[30,60,188,250],'Units','normalized',...
                'data', levels,...
                'ColumnEditable', [false(1,1),true(n_condition_name,1)'],...
                'columnFormat', columnformat,...
                'rowName',[],'RowStriping','off',...
                'ColumnName',columnname,'ColumnWidth','auto');
        else
            reset_levels = levels;
            reset_levels(:,2) = num2cell(false(n_condition_name,1));
            fixed_predictors_table = uitable('position',[30,60,188,250],'Units','normalized',...
                'data', reset_levels,...
                'ColumnEditable', [false(1,1),true(n_condition_name,1)'],...
                'columnFormat', columnformat,...
                'rowName',[],'RowStriping','off',...
                'ColumnName',columnname,'ColumnWidth','auto');
        end
    end

    function option_model(~,~)
        prompt={'singlepredi', 'parallelname', 'LinearMixed Model Optional. For more information please check http://www.mathworks.com/help/stats/fitlme.html#namevaluepairarguments'};
        name='Model Option';
        numlines = [1 15;1 15;5 100];
        text =strcat('''DummyVarCoding''',',','''effect''',',','''FitMethod''',',','''REML''');
        defaultanswer={'1','',text,'type help fitlme for more information'};
        options=inputdlg(prompt,name,numlines,defaultanswer);
        if isempty(options)==0
        opt.singlepredi = str2double(options{1});
        if ~isempty(options{2})
            opt.parallelname = options{2};
        end
        tmp=strsplit(options{3},''',''');
        tmptx=tmp{1};tmp{1}=tmptx(2:end);tmptx=tmp{end};tmp{end}=tmptx(1:end-1);
        opttext=tmp;
        end
    end

    function see_formula(~,~)
        %uiresume(gcf)
        retrieve_random_table = get(random_predictors_table,'data');
        retrieve_fixed_table = get(fixed_predictors_table,'data');
        random_factors = retrieve_random_table(cell2mat(retrieve_random_table(:,2))==1,1);
        fixed_factors = retrieve_fixed_table(cell2mat(retrieve_fixed_table(:,2))==1,1);
        
        random_part = cell(1,length(random_factors));
        for j1=1:length(random_factors)
%             if j1==length(random_factors)
                random_part{1,j1} = strcat(random_factors{j1},')');
%             else
%                 random_part{1,j1} = strcat(random_factors{j1},')+(1|');
%             end
        end
        random_part = strcat('+(1|',random_part);
        fixed_part = implode(fixed_factors,'+');
        formulanew = strcat('PixelIntensity ~ ',' ',fixed_part,strjoin(random_part));
        prompt={''};
        name='Please click the formula before entering into the model';
        numlines = [5 100];
        defaultanswer={formulanew};
        formula=inputdlg(prompt,name,numlines,defaultanswer);
        if isempty(formula)==0
        formula = formula{1};
        % h= msgbox(formula);
        % set(h,'color','white');
        idx_formula = 1; %user edit the formula
        end
    end
end