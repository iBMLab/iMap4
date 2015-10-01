function [answer,Mask,fixthres] = select_mask_parameter(FixMap,Mask,handles)
answer = [];
fixthres=handles.fixthres;
fixthres = round(fixthres*100)/100;
size_figure = [600,300];
S.fh = figure('menubar','none','resize','off','NumberTitle','off','Name','mean fixation bias and the correspondent mask','position',[300 200 size_figure(1) size_figure(2)]);
subplot(1,2,1)
imagesc(squeeze(mean(FixMap,1)));axis('equal','off')
subplot(1,2,2)
imagesc(Mask);
axis('equal','off')

set(gcf,'color','w'); % white background
scale_string = {'Enter a new threshold to regenerate Mask or press done to continue'};
S.text_1 = uicontrol('style','text',...
    'units','pix',...
    'position',[10 70 590 30],...
    'string',scale_string  ,...imagesc(squeeze(mean(handles.FixMap_estimated,1)));axis('equal','off')
    'HorizontalAlignment','center','Fontweight','Bold','Background','White');
S.text = uicontrol('style','text',...
    'units','pix',...
    'position',[130 50 100 20],...
    'string',' Threshold',...
    'HorizontalAlignment','center');

S.edit = uicontrol('style','edit',...
    'unit','pix','string',num2str(fixthres),'backgroundcolor','white','position',[240 50 60 20]);

% create push button to validate the mask threshold value
sz_pshb1 = [60,20];
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Validate','position',[310 50 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',@validate);

% create push button to validate the mask threshold value
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Done','position',[310 20 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',@done);

uiwait(gcf)

% local callback functions
    function validate(hObject,eventdata)
        
        uiresume(gcf)
        fixthres = str2double(get(S.edit,'string'));
        Mask = squeeze(mean(FixMap,1))>fixthres;
        subplot(1,2,1)
        imagesc(squeeze(mean(FixMap,1)));axis('equal','off')
        subplot(1,2,2)
        imagesc(Mask);
        axis('equal','off')
        
        uiwait(gcf)
    end
%---------------------------------------------------------------
    function done(source,event)
        % get values of checkboxes and extract selected predictors
        uiresume(gcf)
        answer = 0;
        if ishandle(S.fh)
            delete(S.fh)
        end
    end
end





