function [himgsc,hStrings] = imsqrmat(mat, varargin)
% pretty display for small (near)square matrix
%--------------------------------------------------------------------------
% Copyright (C) iMap Team 2016

if nargin == 1
    TickX = [];
    TickY = [];
elseif nargin == 2
    TickX = varargin{1};
    TickY = varargin{1};
else
    TickX = varargin{1};
    TickY = varargin{2};
end

himgsc      = imagesc(mat);
axis square;
% display result
textStrings = num2str(mat(:),'%0.2f');  % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x,y]       = meshgrid(1:size(mat,2),1:size(mat,1));   % Create x and y coordinates for the strings
hStrings    = text(x(~isnan(mat(:))),y(~isnan(mat(:))),textStrings(~isnan(mat(:))),...      % Plot the strings
    'HorizontalAlignment','center');
midValue    = mean(get(gca,'CLim'));  % Get the middle value of the color range
textColors  = repmat(mat(~isnan(mat(:))) < midValue,1,3);  % Choose white or black for the
%   text color of the strings so
%   they can be easily seen over
%   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  % Change the text colors
set(gca,'XTick',  1:size(mat,2),...                         % Change the axes tick marks
    'XTickLabel', TickX,...  %   and tick labels
    'YTick',      1:size(mat,1),...
    'YTickLabel', TickY,...
    'TickLength', [0 0],...
    'xdir','reverse','ydir','normal');
end