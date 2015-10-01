function [Predictors_continuous_default,Predictors_categorical_default,Predictors_continuous,Predictors_categorical,counter,handles] = select_conditions(hObject,handles,counter)


% show the original data create table  of the original data

sz_table = [820,400]; %size of table
txtboxheight=30;
%size of label in GUI %sze_label = [120,40];
size_panel = getpixelposition(handles.data_imported,true);
size_panel_button = getpixelposition(handles.panel_button,true);

%create table
t =  uitable(handles.figure1,'position',[size_panel(1)+20,size_panel_button(2)-450,sz_table(1),sz_table(2)],'Units','normalized',...
    'Data', handles.data,'ColumnName', handles.colnames,'RowName',[],'ColumnWidth','auto',...
    'tag','table','RowStriping','off');

% tableextent = get(t,'Extent');
% oldposition = get(t,'Position');
% if tableextent(3)<0.75
%     newposition = [oldposition(1) oldposition(2) tableextent(3) oldposition(4)];
%     set(t, 'Position', newposition);
% end

handles.table1 = t;

% uitable(handles.figure1,'position',[120,40,750,500],'Units','normalized',...
%      'Data', handles.data,'ColumnName',handles.colnames,'RowName',[],'ColumnWidth','auto', 'RowStriping','off');

% Initialisation of some parameters
Predictors_continuous_default = struct; Predictors_categorical_default = struct;
Predictors_continuous =struct;Predictors_categorical=struct;

columns = handles.colnames;
%columns = handles.colnames;
data =  handles.data; % always use the original data.
% create main window
size_figure = [650,500];
S.fh = figure('units','pixels',...
    'position',[300 200 size_figure(1) size_figure(2)],...
    'menubar','none',...
    'name','Predictors',...
    'numbertitle','off',...
    'resize','off');

% movegui(S.fh,'center')
set(gcf,'color','w'); % white background



% create panel to determine the name of columns
sz_panel  = [550,100];
% xpos_panel = ceil((size_figure(1)-sz_panel(2))/2); % center the figure on the screen horizontally
uipanel('parent',S.fh, 'units','pixels', 'pos',[40 350 sz_panel(1) sz_panel(2)]);

% Determine the Subject the Trials X Y and the duration
r = {'Please select the column name of the following predicotrs:'};
S.text_general = uicontrol('style','text','string',r,'position',[20 460 630 txtboxheight],...
    'HorizontalAlignment','center','FontWeight','Bold','Backgroundcolor','white');

rr = {'Subject','Trial','Eye_X','Eye_Y','Duration'};
n_column = {'1','2','3','4','7'}; % by default the removed column

for i = 1:length(n_column)
    
    S.text(i) = uicontrol('style','text',...
        'units','pix',...
        'position',[90*i 420 70 20],...
        'string',rr(i),...
        'HorizontalAlignment','center');
    S.edit(i) = uicontrol('style','popupmenu',...
        'unit','pix',...
        'string',columns,...
        'backgroundcolor','white',...
        'position',[90*i 390 80 20]);
    %'string',num2str(n_column{i}),...
    %
end


% Create push button 'validate_column'

sz_pshb1 = [70,20];
%xpos_pshb1 = ceil((sz_panel(1)-sz_pshb1(2))/2); % center the figure on the screen horizontally
%ypos_table = ceil((size_figure(2)-sz_table(1))/2); % center the figure on the screen vertically

pshb1 = uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Validate','position',[270 360 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',@validate_column);
uiwait(gcf)

if ishandle(S.fh)==0
    % User clicked cancel.
    return;
end





n_column = cell2mat(get(S.edit,'Value')); % update the columns selected by the user
predictors_categorical_default = columns(n_column(1:2))'; %default categorical conditions
idx_categorical_default = n_column(1:2);
predictors_continuous_default  = columns(n_column(3:5))'; %default continuous conditions
idx_continuous_default = n_column(3:5);
Predictors_categorical_default = dataset(predictors_categorical_default,idx_categorical_default);
Predictors_continuous_default = dataset(predictors_continuous_default,idx_continuous_default);

%Determine the Subject the Trials X Y and the duration
r = {'Please double check if the predictor types are correct (continuous or categorical):'};
S.text_general = uicontrol('style','text','string',r,'position',[20 360-75 630 txtboxheight],...
    'HorizontalAlignment','center','FontWeight','Bold','Backgroundcolor','white');

% Column names and column format
columnname = {'   Predictors   ','   CheckBox   ','  Predictors'' type   '};
columnformat = {'char','logical',{'Continuous' 'Categorical'}};

% Define the data
new_columns = columns'; new_columns(n_column) = [];
cell_continuous = cell(length(new_columns),1);

cell_continuous(:) = cellstr('Continuous');

%Determine if a column is categorical or not!
for h=1:length(new_columns)
    idx_columns = strcmp(new_columns{h},handles.colnames)==1;
    vector = unique(handles.data(:,idx_columns)); % to be chekced if we keep original or not
    if length(vector)<=4
        idx_columns_ = strcmp(new_columns{h},new_columns)==1;
        cell_continuous(idx_columns_) =  cellstr('Categorical');
    end
end

d = [new_columns,num2cell(false(length(new_columns),1)), cell_continuous];
clr = '#F0F0F0';
d_vis = d;
d_vis(:,1) = strcat(...
    ['<html><body bgcolor="' clr '" text="#000000" width="100px">'],d_vis(:,1));

% Create the uitable
condition_table = uitable('Data', d_vis,...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', [false true true],...
    'RowName',[],'Position',[180 70 285 200],...
    'Units','normalized',...
    'RowStriping','off');


counter_categorical =  1;    counter_continuous  =  1;

while(counter_categorical==1 && counter_continuous  ==  1)
    if ishandle(S.fh)==1
        % create push button 'validate'
        uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
            'String','Continue','position',[270 40 70 20],...
            'Callback',@validate);
        
        uiwait(gcf)
    else
        counter =0; %user close the figure
        return;
    end
end



%---------------------------------------------------------------
% local callback functions
    function   validate_column(source,event)
        % get values of checkboxes and extract selected predictors
        uiresume(gcf)
        set(pshb1,'enable','off');
        %set(findall(panel1, '-property', 'enable'), 'enable', 'off')
        
    end
% local callback functions
    function validate(source,event)
        % get values of checkboxes and extract selected conditions as well as
        % their type
        uiresume(gcf)
        counter = 1;
        %retrieve all the table
        retrieve_data = get(condition_table,'data');
        
        %retrieve the values of checkbox to detect the choice of the user
        selected_conditions = find(cell2mat(retrieve_data(:,2))==1);
        handles.selected_conditions = length(selected_conditions);
        
        
        
        if isempty(selected_conditions)==0
            % get the predictor type
            
            for n_condition = 1:length(selected_conditions)
                p_type = retrieve_data{selected_conditions(n_condition),3};
                
                switch p_type
                    
                    case 'Categorical'
                        
                        predictors_categorical(counter_categorical,1) = d(selected_conditions(n_condition),1);
                        idx_categorical(counter_categorical,1) = find(strcmp(predictors_categorical(counter_categorical,1),columns)==1); %retrieve the column number of the categorical conditions
                        %predictors_categorical(counter_categorical,2) = idx;
                        counter_categorical = counter_categorical + 1;
                        
                    case 'Continuous'
                        
                        predictors_continuous(counter_continuous,1) = d(selected_conditions(n_condition),1);
                        idx_continuous(counter_continuous,1) = find(strcmp(predictors_continuous(counter_continuous,1),columns)==1); %retrieve the column number of the continuous conditions
                        %predictors_continuous(counter_continuous,2) = idx_continuous;
                        counter_continuous = counter_continuous + 1;
                end
            end
            
            if counter_continuous==1 % the user did not select any continuous conditions Predictors_continuous = dataset(predictors_continuous,idx_continuous);
                
                predictors_continuous = [];
                idx_continuous = [];
                Predictors_continuous = dataset(predictors_continuous,idx_continuous);
                
            else
                Predictors_continuous = dataset(predictors_continuous,idx_continuous);
            end
            
            switch counter_categorical
                
                case 1  %the user did not select any categorical conditions
                    
                    predictors_categorical = [];
                    idx_categorical = [];
                    Predictors_categorical = dataset(predictors_categorical,idx_categorical);
                    
                case 2  % the user selects only one categorical condition
                    
                    Predictors_categorical = dataset(predictors_categorical,idx_categorical);
                    categorical_column = unique(data(:,Predictors_categorical.idx_categorical(1)));  %only one categorical selected, verify if it has at least 2 levels
                    
                    if length(categorical_column)==1
                        counter_categorical =  1;    counter_continuous  =  1;
                        errordlg('You need to select a categorical predictor with at least 2 levels to proceed, or can you can select at least two categorical predictors');
                    end
                    
                otherwise % the user selects more than one categorical condition
                    Predictors_categorical = dataset(predictors_categorical,idx_categorical);
                    %    answer = check_level_categorical_conditions(Predictors_categorical);
                    %    t = 1;
                    
            end
            
            
        else
            
            %the user did not select neither continuous nor categorical
            errordlg('You need to select at least one continuous or categorical condition to proceed');
        end
        
    end
delete(gcf);
end






