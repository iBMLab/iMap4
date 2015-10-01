function [answer,new_data] = check_columns_file(handles)

% Define the data
colnames = handles.colnames_original;
colnames(end) =cellstr('Name of file(s)');

if handles.subject_in ==1
    colnames(end-1) =cellstr('Subject');
    
end
answer = handles.colnames;
new_data = handles.data_original;
d = [handles.colnames_original' colnames'];

clr = '#F0F0F0';
d_vis = d;
d_vis(:,1) = strcat(...
    ['<html><body bgcolor="' clr '" text="#000000" width="100px">'],d_vis(:,1));
size_figure = [370,370];
S.fh = figure('units','pixels',...
    'position',[300 200 size_figure(1) size_figure(2)],...
    'menubar','none',...
    'name','Please select your data file(s)',...
    'numbertitle','off',...
    'resize','off');
movegui(S.fh,'center')
set(gcf,'color','w');

r = {'Please select the column name files','of the imported data :'};
S.text_general = uicontrol('parent',S.fh,'style','text','string',r,'position',[0 330  size_figure(2) 30],...
    'HorizontalAlignment','center','FontWeight','Bold','Backgroundcolor','white');


% size and position of the table
sz_table = [228,280];
% sz_table = [300,380];

% Column names and column format
columnname = {'   Column Index       ','   Column Name       '};
columnformat = {'char','char'};
column_table = uitable('Data', d_vis,...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', [false true],...
    'RowName',[],'Position',[80 45 sz_table(1) sz_table(2)],...
    'Units','normalized',...
    'RowStriping','off');

%pushbutton validate
sz_pshb1 = [70,20];
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Import Column Name','position',[180 15 120 20],...
    'Callback',@import_headers);

uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Continue','position',[100 15 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',@validate_column);


uiwait(gcf)
    function import_headers(~,~)
        uiresume(gcf)
        [filename, pathname] = uigetfile('*.txt',' Select the file that include headers');
        if ischar(filename)
            filetxt=textscan(fopen(strcat(pathname,filename)),'%s');
            headers = filetxt{1};
            
            if handles.subject_in ==1
                if length(headers)~= length(d_vis) + 2
                    h = errordlg('Please verify that the number of column header is correct');
                    set(h,'color','white');
                    uiwait(gcf)
                else
                    d_vis(1:length(headers)-2,2) = headers;
                end
            else
                if length(filetxt{1})~= length(d_vis) + 1
                    d_vis(1:length(d_vis)-1,2) = headers;
                else
                    h = errordlg('Please verify that the number of column header is correct');
                    set(h,'color','white');
                    uiwait(gcf)
                end
            end
            
            column_table = uitable('Data', d_vis,...
                'ColumnName', columnname,...
                'ColumnFormat', columnformat,...
                'ColumnEditable', [false true],...
                'RowName',[],'Position',[80 45 sz_table(1) sz_table(2)],...
                'Units','normalized',...
                'RowStriping','off');
            uiwait(gcf)
        end
    end

    function   validate_column(~,~)
        % get name of columns and retrieve all the table
        uiresume(gcf)
        
        retrieve_data = get(column_table,'data');
        answer  = retrieve_data(:,2)';
        if ishandle(S.fh)
            delete(S.fh)
        end
        
        %Ask the user if he wants to combine categorical conditions
        answer_user=questdlg('Would you like to create new predictors?','Create New Predictors','yes','no','no');
        if strcmp(answer_user,'yes')
            %Determine the Subject the Trials X Y and the duration
            size_figure = [370,370];
            S.fh = figure('units','pixels',...
                'position',[300 200 size_figure(1) size_figure(2)],...
                'menubar','none',...
                'name','Predictors',...
                'numbertitle','off',...
                'resize','off');
            movegui(S.fh,'center')
            set(gcf,'color','w');
            r = {'Please select the predictors you want to combine'};
            S.text_general = uicontrol('style','text','string',r,'position',[0 330  size_figure(2) 30],...
                'HorizontalAlignment','center','FontWeight','Bold','Backgroundcolor','white');
            
            % Create uitable to determine the name of columns
            %n=length(columns)-length(n_column); % the number of the rest of columns
            
            % Column names and column format
            columnname = {'   Predictors   ','   CheckBox   '};
            columnformat = {'char','logical'};
            
            % Define the data
            new_columns = answer';
            d = [new_columns num2cell(false(length(new_columns),1))];
            
            clr = '#F0F0F0';
            d_vis = d;
            d_vis(:,1) = strcat(...
                ['<html><body bgcolor="' clr '" text="#000000" width="100px">'],d_vis(:,1));
            
            % size and position of the table
            sz_table = [178,200];
            xpos_table = ceil((size_figure(1)-sz_table(2))/2); % center the figure on the screen horizontally
            
            % Create the uitable
            predictors_table =  uitable('Data', d_vis,...
                'ColumnName', columnname,...
                'ColumnFormat', columnformat,...
                'ColumnEditable', [false true],...
                'RowName',[],'Position',[xpos_table 120 sz_table(1) sz_table(2)],...
                'Units','normalized',...
                'RowStriping','off');
            
            % Create push button to 'combine predictors'
            uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
                'position',[80 50 100 20],...
                'String', 'Combine Predictors', 'Callback',{@combine_condition,handles});
            
            % Create push button to 'cancel'
            uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
                'position',[200 50 70 20],...
                'String', 'Cancel', 'Callback',@cancel_combine);
            uiwait(gcf)
        end
        if ishandle(S.fh)
            delete(S.fh)
        end
        
        function combine_condition(~,~,handles)
            
            uiresume(gcf)
            %retrieve all the table
            retrieve_data = get(predictors_table,'data');
            %retrieve the values of checkbox to detect the choice of the user
            selected_conditions = find(cell2mat(retrieve_data(:,2))==1);
            handles.n_newconditions = selected_conditions; %selected columns to combine
            handles.selected_conditions = length(selected_conditions);
            condition_not_treated  = 0;
            for j = 1:length(selected_conditions)
                
                %Display table for user
                categorical_column  = find(ismember(new_columns,d{selected_conditions(j),1})); % find the categorical column
                Levels = unique(handles.data(:,categorical_column)); % retrieve the levels of the categorical columns
                
                if length(Levels)<=2
                    condition_not_treated = condition_not_treated+1; % if the condition has few levels (less or equal to 2),the user cannot combine levels
                    condition_name = d{selected_conditions(j),1};
                    msg = msgbox(strcat('The following predictor',{' '},'"',condition_name,'"',{' '},'has only two levels or less, for that reason it is not possible to combine levels.'));
                    set(msg, 'color', 'white');
                    uiwait(msg)
                    continue;
                else
                    
                    
                    levels = [Levels,num2cell(false(length(Levels),length(Levels)-1))];
                    n_levels = length(levels(:,1))-1;
                    
                    
                    %set the column according to the number of levels
                    column_name = char2cell(num2str(1:n_levels))';
                    columnname_categorical = [d{selected_conditions(j),1},column_name];
                    column_type = cell(n_levels,1);
                    column_type(:) = cellstr('logical');
                    columnformat_categorical = ['char',column_type'];
                    screensize = get(0,'ScreenSize');
                    %sz = [185,180];
                    sz = [500,300];
                    xpos = ceil((screensize(3)-sz(2))/2); % center the figure on the screen horizontally
                    ypos = ceil((screensize(4)-sz(1))/2); % center the figure on the screen vertically
                    
                    f = figure('Name','Combine Levels','Position',[xpos ypos sz(1) sz(2)],...
                        'menubar','none','numbertitle','off','resize','off','units','normalized');
                    set(gcf,'color','w'); % white background
                    
                    % Create the uitable
                    categorical_table = uitable('parent',f,'position',[0,50,500,250],'Units','normalized',...
                        'data', levels,...
                        'ColumnEditable', [false(1,1),true(n_levels,1)'],...
                        'columnFormat', columnformat_categorical,...
                        'rowName',[],'RowStriping','off',...
                        'ColumnName',columnname_categorical);
                    
                    % create push button 'validate'
                    uicontrol('parent',f,'Style','pushbutton','Units','pixel',...
                        'String','Validate','position',[150 10 70 20],...
                        'Callback', {@validate_levels, categorical_table,categorical_column});
                    
                    % create push button 'cancel'
                    uicontrol('parent',f,'Style','pushbutton','Units','pixel',...
                        'String','Cancel','position',[240 10 70 20],...
                        'Callback',@cancel_levels);
                    uiwait(gcf)
                end
            end
            %---------------------------------------------------------------
            % local callback functions
            
            function validate_levels(~,~,categorical_table,categorical_column)
                
                %retrieve all the table
                retrieve_categorical_table = get(categorical_table,'data');
                prompt = {'Enter the name of this new predictor'}; % ask the user to name the level condition
                dlg_title = ''; num_lines= 1; def     = {''};
                new_level_name  = inputdlg(prompt,dlg_title,num_lines,def); % the name of the level
                
                
                %retrieve the values of checkbox to detect the levels selected by the user
                all_selected_levels = cell2mat(retrieve_categorical_table(:,2:end));
                [row,all_col] = find(all_selected_levels==1);
                col = unique(all_col);
                last_column_data = size(handles.data,2);
                
                new_condition = handles.data(:,categorical_column);
                
                for k = 1:length(col)
                    selected_levels = cell2mat(retrieve_categorical_table(:,col(k)+1))==1; % we add + 1 because we do not consider the first column which is the name of the parameters
                    level_name = retrieve_categorical_table(selected_levels,1);
                    position_condition =  ismember(new_condition,level_name)==1;
                    if length(level_name)>1
                        new_condition(position_condition) = cellstr(strjoin(level_name','&')); % replace levels of conditions with the one named by the user
                    else
                        
                        new_condition(position_condition) = cellstr(level_name); % replace levels of conditions with the one named by the user
                        
                    end
                    
                    
                end
                
                handles.data(:,last_column_data+1) = new_condition; % replace levels of conditions with the one named by the user
                
                %remove the non-selected parameters
                IDX = [];
                all_non_selected_levels = setdiff(1:length(retrieve_categorical_table(:,1)),row);
                for k=1:length(all_non_selected_levels)
                    idx = find(strcmp(retrieve_categorical_table(all_non_selected_levels(k),1),handles.data(:,last_column_data+1))==1);
                    IDX = [IDX;idx];
                end
                handles.data(IDX,last_column_data+1) = cellstr('');
                answer(last_column_data+1) = new_level_name(1);
                handles.colnames = answer;
                delete(gcf)
                new_data = handles.data;
                
                sz_table = [820,400]; %size of table
                % size of label in GUI %sze_label = [120,40];
                size_panel = getpixelposition(handles.data_imported,true);
                size_panel_button = getpixelposition(handles.panel_button,true);
                uitable(handles.figure1,'position',[size_panel(1)+20,size_panel_button(2)-450,sz_table(1),sz_table(2)],...
                    'Data', new_data,'ColumnName', answer,'RowName',[],'ColumnWidth','auto',...
                    'Units','normalized','tag','table','RowStriping','off');
                
            end
            
            function cancel_levels(~,~)
                
                delete(gcf)
            end
            
            
        end
    end
    function cancel_combine(~,~)
        delete(gcf)
    end
if ishandle(S.fh)
    delete(S.fh)
end
end