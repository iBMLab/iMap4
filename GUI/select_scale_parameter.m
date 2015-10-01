
function [scale_parameter, ySize2, xSize2,idx] = select_scale_parameter(xSize,ySize)
idx = 0;
scale_parameter = round(25000/max(xSize,ySize))/100;
[ySize2,xSize2]=size(imresize(ones(ySize,xSize),scale_parameter)); % to verify!
scale_parameter = num2str(scale_parameter);

screensize = get(0,'ScreenSize');
size_figure = [350,160];

xpos = ceil((screensize(3)-size_figure(2))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-size_figure(1))/2); % center the figure on the screen vertically


S.fh = figure('units','pixels',...
    'position',[xpos ypos size_figure(1) size_figure(2)],...
    'menubar','none',...
    'name','Scale Parameter',...
    'numbertitle','off',...
    'resize','off');

set(gcf,'color','w'); % white background
scale_string = {'Enter the scale parameter and press enter to get' 'the updated values of X and Y size'};
S.text_1 = uicontrol('style','text',...
    'units','pix',...
    'position',[10 125 340 30],...
    'string',scale_string  ,...
    'HorizontalAlignment','center','Fontweight','Bold','Background','White');
xsize_value = strcat('  Original  X size =', {' '},num2str(xSize));
xsize_string = {'  Updated X size' xsize_value{1}};

S.text_2 = uicontrol('style','text',...
    'units','pix',...
    'position',[40 95 130 20],...
    'string',' Scale parameter',...
    'HorizontalAlignment','left');

S.text_3 = uicontrol('style','text',...
    'units','pix',...
    'position',[40 55 130 30],...
    'string',xsize_string,...
    'HorizontalAlignment','left');

ysize_value = strcat('  Original  Y size =', {' '},num2str(ySize));
ysize_string = {'  Updated Y size' ysize_value{1}};
S.text_4 = uicontrol('style','text',...
    'units','pix',...
    'position',[40 15 130 30],...
    'string',ysize_string,...
    'HorizontalAlignment','left');

S.edit_2 = uicontrol('style','edit',...
    'unit','pix',...
    'string',xSize2,...
    'backgroundcolor','white',...
    'position',[190 59 80 20],'callback',@edit_callback_2);


S.edit_3 = uicontrol('style','edit',...
    'unit','pix',...
    'string',ySize2,...
    'backgroundcolor','white',...
    'position',[190 19 80 20],'callback', @edit_callback_3);
% Create push button 'validate_column'
sz_pshb1 = [70,20];
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Validate','position',[275 95 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',@validate);

%set(S.fh,'windowbuttonmotionfcn',edit_callback);
S.edit_1 = uicontrol('style','edit',...
    'unit','pix',...
    'string',num2str(scale_parameter),...
    'backgroundcolor','white',...
    'position',[190 95 80 20],'callback', @edit_callback);




uiwait(gcf)


    function edit_callback_2(hObject,~)
        
        
        index_selected = get(hObject,'string');
        scale_parameter = round((str2double(index_selected)/xSize)*100)/100;
        [ySize2,xSize2]=size(imresize(ones(ySize,xSize),scale_parameter));
        set(S.edit_1,'string',scale_parameter)
        %set(S.edit_2,'string',xSize2)
        set(S.edit_3,'string',ySize2)
        
    end

    function edit_callback_3(hObject,~)
        
        
        index_selected = get(hObject,'string');
        scale_parameter = round((str2double(index_selected)/ySize)*100)/100;
        [ySize2,xSize2]=size(imresize(ones(ySize,xSize),scale_parameter));
        set(S.edit_1,'string',scale_parameter)
        set(S.edit_2,'string',xSize2)
        %set(S.edit_3,'string',ySize2)
        
    end

    function edit_callback(hObject,~)
        
        
        index_selected = get(hObject,'string');
        
        [ySize2,xSize2]=size(imresize(ones(ySize,xSize),str2double(index_selected)));
        set(S.edit_2,'string',xSize2)
        set(S.edit_3,'string',ySize2)
        
        % 5    disp(['Called by object with handle: ' num2str(callingObj) ])
        %     disp(['Character pressed: ' data.Character])
        %     disp(['Modifier: ' data.Modifier])
        %     disp(['Key: ' data.Key])
        %     disp(['String in editbox: ' get(callingObj,'string') ] )
        %set(S,'string','a')
        
    end
%---------------------------------------------------------------
% local callback functions
    function   validate(~,~)
        % get values of checkboxes and extract selected predictors
        uiresume(gcf)
        idx = 1;
        scale_parameter = get(S.edit_1,'String');
       delete(S.fh)
    end





end
