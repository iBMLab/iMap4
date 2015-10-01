function answer = select_parameters()

%prompt={'Enter screen x in pixel','Enter screen y in pixel', 'Enter image x in pixel', 'Enter image y in pixel'};

prompt={'Screen x in pixel','Screen y in pixel', 'Image x in pixel',...
    'Image y in pixel','Horizontal screen size in cm','Vertical screen size in cm','Distance between the subject and the screen in cm'};

name= ''; %'Input Parameters';
numlines=1;
scrsz=get(0,'ScreenSize');
x_screen = num2str(scrsz(3));
y_screen = num2str(scrsz(4));
defaultanswer = {x_screen,y_screen,x_screen,y_screen,'','',''};
answer=inputdlg(prompt,name,numlines,defaultanswer);
if isempty(answer)==0
%check of even number of x image and y image
if mod(str2double(answer(3)),2)~=0
   answer(3) = cellstr(num2str(str2double(answer(3))+1)); 
end

if mod(str2double(answer(4)),2)~=0
   answer(4) = cellstr(num2str(str2double(answer(4))+1)); 
end
end

