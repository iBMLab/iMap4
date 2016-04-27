function [answer,Mask,fixthres] = select_mask_parameter(FixMap,Mask,handles)
answer = [];
fixthres=handles.fixthres;
 % This line is removed because if we have a small threshold 0.0025 for
 % example by rounding the number the threshold will be considered as 0!
% fixthres = round(fixthres*100)/100;
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
    'position',[10 20 590 30],...
    'string',scale_string  ,...imagesc(squeeze(mean(handles.FixMap_estimated,1)));axis('equal','off')
    'HorizontalAlignment','center','Fontweight','Bold','Background','White');
S.text = uicontrol('style','text',...
    'units','pix',...
    'position',[130 10 100 20],...
    'string',' Threshold',...
    'HorizontalAlignment','center');

S.edit = uicontrol('style','edit',...
    'unit','pix','string',num2str(fixthres),'backgroundcolor','white','position',[240 10 60 20]);

% create push button to validate the mask threshold value
sz_pshb1 = [60,20];
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Validate','position',[310 10 sz_pshb1(1) sz_pshb1(2)],...
    'Callback',@validate);

% create push button to validate the mask threshold value
uicontrol('parent',S.fh,'Style','pushbutton','Units','pixel',...
    'String','Done','position',[390 10 sz_pshb1(1) sz_pshb1(2)],...
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





