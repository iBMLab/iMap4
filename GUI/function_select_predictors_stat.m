function [choice_conditions,selection_error] = function_select_predictors_stat(ConditionM)
choice_conditions = [];
selection_error = 0;
screensize = get(0,'ScreenSize');
sz = [185,180];
size_figure = [300,300];

%sz = [300,200];
xpos = ceil((screensize(3)-sz(2))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-sz(1))/2); % center the figure on the screen vertically

S.fh = figure('units','pixels','position',[xpos, ypos size_figure(1) size_figure(2)],...
    'menubar','none',...
    'name','Joint Conditions',...
    'numbertitle','off',...
    'resize','off');
%movegui(S.fh,'center')
set(gcf,'color','w');

r = {'Kindly select the categorical conditions to see their corresponding boxplots'};
S.text_general = uicontrol('parent',S.fh,'style','text','string',r,'position',[50 250  200 30],...
    'HorizontalAlignment','center','FontWeight','Bold','Backgroundcolor','white');
if isa(ConditionM,'dataset')
    CName = ConditionM.Properties.VarNames;
elseif isa(ConditionM,'table')
    CName   = ConditionM.Properties.VariableNames;
end
% remove subject from CName
CName(strcmp(CName,'subject')==1)=[];
n = length(CName);
levels = [CName', num2cell(true(n,1))];
%sz_table_random = [160,150];
columnname = {'Categorical|Predictors','CheckBox'};
column_type = cell(1,1);
column_type(:) = cellstr('logical');
columnformat = ['char',column_type'];



categorical_predictor_table = uitable(S.fh,'position',[30,30,226,250],'Units','normalized',...
    'data', levels,...
    'ColumnEditable', [false,true],...
    'columnFormat', columnformat,...
    'rowName',[],'RowStriping','off',...
    'ColumnName',columnname);

sz_pshb1 = [70,20];

uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Continue','position',[40 8 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',{@validate_predictors,categorical_predictor_table});

% create push button 'validate'
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Skip and continue to LMM','position',[120 8 130 20],...
    'Callback',@cancel_predictors);

uiwait(gcf)
if ishandle(S.fh)
    delete(S.fh)
end

    function validate_predictors(~,~,categorical_predictor_table)
        uiresume(gcf)
        retrieve_data = get(categorical_predictor_table,'data');
        choice_conditions = retrieve_data(cell2mat(retrieve_data(:,2))==1,1);
        if length(choice_conditions)<2
            selection_error = 1;
            h1 = errordlg('Please select at least two predictors to proceed');
            uiwait(h1);
        end
    end
    function cancel_predictors(~,~)
        delete(gcf)
        choice_conditions = [];
    end
end