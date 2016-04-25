function [data,counter,mapType,check_number_fixations] = function_find_fixation(handles,fix_number)
data = handles.data; %get the data
counter = 0;
mapType = 1;
check_number_fixations = 1; % this indicates that we have more than one fixation per trial and thus would able to do the estimate and the single trial method

% find the fixation index and and the maximumn number of fixations
% retrieve the trial column data
categorical_conditions = handles.predictors_categorical_default;


%Transform the trial into numbers
 [C, ic,trial] = unique(handles.data(:,categorical_conditions.idx_categorical_default(2)),'stable');
if fix_number~=1
if(any(isnan(trial)==1)) % if this condition is true, something is wrong with the trial column
    max_number_fixation = 1;
    fixation_index = [];
else
    % find the number of fixation
    number_fixation = findseq(trial);
    
    % find the maximum number of fixation
    max_number_fixation = max(number_fixation(:,4));
       
    
    % initialize the vector
    fixation_index = zeros(length(number_fixation),1);
    
    % a loop to create the fixation index
    for j = 1:size(number_fixation,1)
        fixation_index(number_fixation(j,2):number_fixation(j,3)) = 1:number_fixation(j,4);
    end
end
else
    counter = 1;
    mapType = 2;
    max_number_fixation = -1;% -1 is just to distiguish between one fixation or when the user enter a wrong column for trial!
    fixation_index = ones(1,length(trial));
end

if isempty(max_number_fixation)==1
    errordlg('Unable to retrieve the number of fixation. Make sure you select the right column for trial (i.e., numerical values).');
    
elseif max_number_fixation>1
    size_figure = [340,150];
    screensize = get(0,'ScreenSize');
    xpos_1 = ceil((screensize(3)-size_figure(2))/2); % center the figure on the screen horizontally
    ypos_1 = ceil((screensize(4)-size_figure(1))/2); % center the figure on the screen vertically
    %%
    %300,200
    S.fh = figure('units','pixels',...
        'position',[xpos_1 ypos_1 size_figure(1) size_figure(2)],...
        'menubar','none',...
        'name','Fixation Options',...
        'numbertitle','off',...
        'resize','off');
    set(S.fh,'color','white')
    %movegui(S.fh,'center');
    % type fixation
    type_fixation0 = {'Duration Map','Number of Fixation Map'};
    % create checkboxes
    hBtnGrp0 = uibuttongroup('Position',[-.01 0 1.02 1],'BackgroundColor', [1 1 1]);
    %set(hBtnGrp,'color','white');
    
    boxhandles0(1) =  uicontrol('Style','Radio', 'Parent',hBtnGrp0, 'HandleVisibility','off',...
        'Position',[50 100  100 20], 'String',type_fixation0{1}) ;
    
    boxhandles0(2) =  uicontrol('Style','Radio', 'Parent',hBtnGrp0, 'HandleVisibility','off',...
        'Position',[160 100  150 20], 'String',type_fixation0{2}) ;
    
    r = {'What type of fixation map would you like?'};
    S.text_general = uicontrol('style','text','string',r,'position',[10 125 300  15],...
        'HorizontalAlignment','center','FontWeight','Bold','Backgroundcolor','white');
    %
    % type fixation
    type_fixation = {'Specific','Exclude','All fixations'};
    % create checkboxes
    hBtnGrp = uibuttongroup('Position',[-.01 -0.01 1.02 .6],'BackgroundColor', [1 1 1]);
    
    %     for k=1:3
    %         boxhandles(k) =  uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off',...
    %             'Position',[30+ 100*(k-1) 30  80 20], 'String',type_fixation{k}) ;
    %     end
    
    boxhandles(1) =  uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off',...
        'Position',[30  50  80 20], 'String',type_fixation{1}) ;
    
    boxhandles(2) =  uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off',...
        'Position',[130 50  80 20], 'String',type_fixation{2}) ;
    
    boxhandles(3) =  uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off',...
        'Position',[230 50  80 20], 'String',type_fixation{3}) ;
    
    r = {'  Do you want to use/exclude specific fixations conditions?'};
    S.text_general = uicontrol('style','text','string',r,'position',[10 75 325 15],...
        'HorizontalAlignment','left','FontWeight','Bold','Backgroundcolor','white');
    
    sz_pshb1 = [70,20];
    uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
        'String','Validate','position',[130 3.5 sz_pshb1(1) sz_pshb1(2)],...
        'Callback',{@validate_type_fixation,boxhandles,boxhandles0,max_number_fixation,fixation_index});
    set(hBtnGrp,'SelectedObject', boxhandles(3));
    uiwait(gcf)
    %     if ishandle(S.fh)==0
    %         % User closed the popup window.
    %         counter = 0;
    %         return;
    %     end
    
else
   check_number_fixations = 0; % This means we have only one fixation and thus only the estimated method should be used 
   size_figure = [340,100];
    screensize = get(0,'ScreenSize');
    xpos_1 = ceil((screensize(3)-size_figure(2))/2); % center the figure on the screen horizontally
    ypos_1 = ceil((screensize(4)-size_figure(1))/2); % center the figure on the screen vertically
    %%
    %300,200
    S.fh = figure('units','pixels',...
        'position',[xpos_1 ypos_1 size_figure(1) size_figure(2)],...
        'menubar','none',...
        'name','Fixation Options',...
        'numbertitle','off',...
        'resize','off');
    set(S.fh,'color','white')
    %movegui(S.fh,'center');
    % type fixation
    type_fixation0 = {'Duration Map','Number of Fixation Map'};
    % create checkboxes
    hBtnGrp0 = uibuttongroup('Position',[-.01 0 1.02 1],'BackgroundColor', [1 1 1]);
    %set(hBtnGrp,'color','white');
    
    boxhandles0(1) =  uicontrol('Style','Radio', 'Parent',hBtnGrp0, 'HandleVisibility','off',...
        'Position',[30 60  100 20], 'String',type_fixation0{1}) ;
    
    boxhandles0(2) =  uicontrol('Style','Radio', 'Parent',hBtnGrp0, 'HandleVisibility','off',...
        'Position',[160 60  150 20], 'String',type_fixation0{2}) ;
    
    r = {'What type of fixation map would you like?'};
    S.text_general = uicontrol('style','text','string',r,'position',[10 125 300  15],...
        'HorizontalAlignment','center','FontWeight','Bold','Backgroundcolor','white');
 
    boxhandles = []; % we go to the case where we use all fixations since we only have one fixation!
    sz_pshb1 = [70,20];
    uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
        'String','Validate','position',[120 20 sz_pshb1(1) sz_pshb1(2)],...
        'Callback',{@validate_type_fixation,boxhandles,boxhandles0,max_number_fixation,fixation_index});

    uiwait(gcf)
    
     
  
  
end

% local callback functions
    function validate_type_fixation(~,~,boxhandles,boxhandles0,max_number_fixation,fixation_index)
        
        % get values of checkboxes and extract selected predictors
        uiresume(gcf)
        
        %retrieve the choice of the uster
        choice_user = find(cell2mat(get(boxhandles,'Value'))==1);
        if isempty(choice_user)==1 % We are in the case of one fixation
        choice_user = 3; % we need all data!
        end
        %retrieve the choice of the uster
        mapType = find(cell2mat(get(boxhandles0,'Value'))==1);

        %delete(S.fh);
        switch choice_user
            
            case {1,2} %user select specific fixations
                
                %sz = [185,180];
                sz = [250,250];
                xpos = ceil((screensize(3)-sz(2))/2); % center the figure on the screen horizontally
                ypos = ceil((screensize(4)-sz(1))/2); % center the figure on the screen vertically
                
                % create figure
                f = figure('Name','Fixations','Position',[xpos ypos sz(1) sz(2)],...
                    'menubar','none','numbertitle','off','resize','off','units','normalized');
                set(gcf,'color','w');
                
                nfixation = max_number_fixation;
                fixations = cell(nfixation,2);
                for i = 1:nfixation
                    fixations(i,:) =[cellstr(strcat('Fixation',{' '}, num2str(i))),num2cell(false,1)'];
                end
                
                columnformat_fixations = {'char','logical'};
                
                % create the uitable
                fixation_table = uitable('parent',f,'position',[10,40,230,200],'Units','normalized',...
                    'data', fixations,...
                    'ColumnEditable', [false,true],...
                    'columnFormat', columnformat_fixations,...
                    'rowName',[],'RowStriping','off',...
                    'ColumnName',{'                   Fixations              ',' Checkbox     '});
                
                % create push button 'validate'
                uicontrol('parent',f,'Style','pushbutton','Units','pixel',...
                    'String','Validate','position',[40 5 70 20],...
                    'Callback', {@validate_fixations,fixation_table,fixation_index,choice_user});
                
                %create push button 'cancel'
                uicontrol('parent',f,'Style','pushbutton','Units','pixel',...
                    'String','Cancel','position',[120 5 70 20],...
                    'Callback', @cancel_fixations);
                uiwait(gcf)
                
            case 3 %user choose all fications (Default)
                %nothing to do here  we use all fixations
                data = handles.data; % get data
                counter = 1;
                if ishandle(S.fh)
                    delete(S.fh);
                end
                
                                
        end
        
        % local callback functions
        function validate_fixations(~,~,fixation_table,fixation_index,choice_user)
            % get values of checkboxes and extract selected predictors
            uiresume(gcf)
            switch  choice_user
                case 1 %use specific fixations
                    all_fixations = get(fixation_table,'data');
                    selected_fixations = find(cell2mat(all_fixations(:,2))==1);
                    if isempty(selected_fixations)==0
                        % find the selected fixations in the data
                        idx_fixation= ismember(fixation_index,selected_fixations);
                        data = data(idx_fixation,:);
                        delete(f)
                        if ishandle(S.fh)
                            delete(S.fh);
                        end
                        counter = 1;
                    else
                        h = errordlg('Error, Pease select at least one fixation');
                        set(h,'color','white');
                    end
                    
                case 2 %use
                    all_fixations = get(fixation_table,'data');
                    selected_fixations = find(cell2mat(all_fixations(:,2))==1);
                    if isempty(selected_fixations)==0
                        % find the selected fixations in the data
                        idx_fixation = ismember(fixation_index,selected_fixations);
                        data(idx_fixation,:) = [];
                        delete(f)
                        if ishandle(S.fh)
                            delete(S.fh);
                        end
                        counter = 1;
                    else
                        h = errordlg('Error, Please select at least one fixation');
                        set(h,'color','white');
                    end
            end
            
        end
        
        function cancel_fixations(~,~)
            delete(gcf)
            uiwait(gcf)
        end
    end

end

