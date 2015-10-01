function [data, subject_in] = create_subject(filename,data,handles)

% prompt={'Enter the number of files that refers to one subject'};
% name= ''; %'Input Parameters';
% numlines=1;
% defaultanswer = {''};
% answer=inputdlg(prompt,name,numlines,defaultanswer);    
size_figure = [350,350];

screensize = get(0,'ScreenSize');
sz = [185,180];
%sz = [300,200];
xpos = ceil((screensize(3)-sz(2))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-sz(1))/2); % center the figure on the screen vertically

S.fh = figure('units','pixels','position',[xpos, ypos size_figure(1) size_figure(2)],...
              'menubar','none',...
              'name','Subject Input',...
              'numbertitle','off',...
              'resize','off');
          %movegui(S.fh,'center')
           set(gcf,'color','w');

r = {'Kindly check the files that belongs to one subject :'};
   S.text_general = uicontrol('parent',S.fh,'style','text','string',r,'position',[0 310  size_figure(2) 30],...
   'HorizontalAlignment','center','FontWeight','Bold','Backgroundcolor','white');
   

%set the column according to the number of levels
if iscell(filename)==0
n_levels = 1;
else
n_levels = length(filename);
end
levels = [filename',num2cell(false(n_levels,n_levels-1))];
column_name = char2cell(num2str(1:n_levels))';
columnname = ['File name',column_name];
column_type = cell(n_levels,1);
column_type(:) = cellstr('logical');
columnformat = ['char',column_type'];



subject_table = uitable('position',[50,60,230,260],'Units','normalized',...
            'data', levels,... 
            'ColumnEditable', [false(1,1),true(n_levels,1)'],...
            'columnFormat', columnformat,...
            'rowName',[],'RowStriping','off',...
            'ColumnName',columnname);

%%%%%%%%%%%%%%%%%%%%%%%%
sz_pshb1 = [70,20];

uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
        'String','Validate','position',[90 15 sz_pshb1(1) sz_pshb1(2)],...
        'Callback',{@validate_levels,subject_table});
    
% create push button 'validate'
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
        'String','Cancel','position',[170 15 70 20],...
        'Callback',@cancel_levels);    
uiwait(gcf)     
    
    
    
    
 function validate_levels(source,event,subject_table)
    
%retrieve all the table
retrieve_subject_table = get(subject_table,'data'); 

all_selected_levels = cell2mat(retrieve_subject_table(:,2:end));
[row,all_col] = find(all_selected_levels==1);
col = unique(all_col);
last_column_data = size(data,2);
new_condition = data(:,last_column_data);

for k = 1:length(col)
selected_levels = cell2mat(retrieve_subject_table(:,col(k)+1))==1; % we add + 1 because we do not consider the first column which is the name of the parameters   
level_name = retrieve_subject_table(selected_levels,1);

%if ischar(level_name{1})
position_condition =  ismember(new_condition,level_name)==1;
if length(level_name)>1
new_condition(position_condition) = cellstr(num2str(k)); % replace levels of conditions with the one named by the user
else

new_condition(position_condition) = cellstr(num2str(k)); % replace levels of conditions with the one named by the user
    
end


end

data(:,last_column_data+1) = new_condition; % replace levels of conditions with the one named by the user
data(:,[end-1,end])=data(:,[end,end-1]);
subject_in = 1;

delete(gcf)

% create table 
% uitable(handles.figure1,'position',[120,40,750,500],'Units','normalized',...
%      'Data', data,'ColumnName', columnames,'RowName',[],'ColumnWidth','auto', 'RowStriping','off');
  end
  
function cancel_levels(source,event)
subject_in = 0;
 delete(gcf)
end

% new_data = [];
% ID_subject = ones(length(filename),1);
% max_loop = ceil(length(filename)/str2double(answer));%know integer for the loop
% jj=0;
% for j=1:max_loop
% ID_subject(jj+1:str2double(answer)+jj) = j;
% jj=jj+str2double(answer);
% end
% for i=1:length(filename)
% load(filename{i});
% variable = fieldnames(load(filename{i}));
% new_data = [new_data;eval(variable{1}),ones(length(eval(variable{1})),1).*ID_subject(i)];
% end
end
