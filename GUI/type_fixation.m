function mapType = type_fixation()
mapType = 1;
size_figure = [340,80];
screensize = get(0,'ScreenSize');
xpos_1 = ceil((screensize(3)-size_figure(2))/2); % center the figure on the screen horizontally
ypos_1 = ceil((screensize(4)-size_figure(1))/2); % center the figure on the screen vertically

%300,200
S.fh = figure('units','pixels',...
    'position',[xpos_1 ypos_1 size_figure(1) size_figure(2)],...
    'menubar','none',...
    'name','Type of Fixation',...
    'numbertitle','off',...
    'resize','off');

% type fixation
type_fixation = {'Duration Map','Number of Fixation Map'};
% create checkboxes
hBtnGrp = uibuttongroup('Position',[0 -0.02  100 400]);
%set(hBtnGrp,'color','white');

boxhandles(1) =  uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off',...
    'Position',[50 30  100 20], 'String',type_fixation{1}) ;

boxhandles(2) =  uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off',...
    'Position',[160 30  150 20], 'String',type_fixation{2}) ;

r = {'What type of fixation map would you like?'};
S.text_general = uicontrol('style','text','string',r,'position',[10 60 300  15],...
    'HorizontalAlignment','center','FontWeight','Bold','Backgroundcolor','white');

sz_pshb1 = [70,20];
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Validate','position',[85 3.5 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',@validate_user_choice);
uiwait(gcf)

    function validate_user_choice(~,~)
        %local callback functions
        
        % get values of checkboxes and extract selected predictors
        uiresume(gcf)
        
        %retrieve the choice of the uster
        mapType = find(cell2mat(get(boxhandles,'Value'))==1);
        close(gcf)
    end
end