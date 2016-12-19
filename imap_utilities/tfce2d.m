function tfce_score = tfce2d(varargin)

% implementation of the Threshold-free cluster enhancement method
% developped for fMRI by Smith & Nichols, NeuroImage 44(2009), 83-98
%
% INPUT tfce_score = tfce2d(data)
%       tfce_score = tfce2d(data,E,H,dh)
%
%       data can be either 2D: a map of t/F values 
%       or data can be 3D: a set of t/F maps computed under H0 
%       E, H and dh are the parameters of the tfce algorithm: defaults are 0.5, 2, 0.1
%       tfce = sum(extend(h)^E*height^H*dh)      
%
% OUPUT tfce_score is a map of scores
%
% ref Pernet, C., Rousselet, G.A., Nichols, T.E. Threshold-free cluster
% enhancement for MEEG: simulation and comparison with the cluster-mass approach.
% NeuroImage, XX, 2012, XX-XX.
%
% Cyril Pernet 18-10-2011
% GAR 21-02-2012: made 2d version for classification images using bwlabel
%                 changed check for negative values to min(data(:)) >= 0  
% modified by Junpeng Lao for iMap4
% -----------------------------
% Copyright (C) LIMO Team 2010


%% default precision for tfce

precision = 100; % define how many thresholds between min t/F map and max t/F map

%% check input

if nargin < 1
    error('not enough input arguments')
elseif nargin == 1 % use default parameters - see Smith & Nichols, NeuroImage 44(2009), 83-98
    E = 0.5; % default for 2D - 0.5 in limo_tfce
    H = 2;
    dh = 0.1;
elseif nargin == 4
    E = varargin{2};
    H = varargin{3};
    dh = varargin{4};
elseif nargin > 5
    error('too many input arguments')
end

data = varargin{1};
[x,y,b]=size(data);
if b == 1
    type = 1;
else
    type = 2;
end

%% start tcfe 

switch type
    
    case{1}
   % ------- tfce real data ----------- 
      
   % define increment size
   increment = (max(data(:)) - min(data(:))) / precision;
   
   % check negative values if so do negate and add scores
   if min(data(:)) >= 0
       
       % select a height, obtain cluster map, obtain extend map (=cluster
       % map but with extend of cluster rather than number of the cluster)
       % then tfce score for that height
       index = 1; tfce = NaN(x,y,length(min(data(:)):increment:max(data(:))));
       for h=min(data(:)):increment:max(data(:))
           [clustered_map, num] = bwlabel(data > h);
           extend_map = zeros(x,y); % same as cluster map but contains extend value instead
           extend_map = integrate(clustered_map,num,extend_map);
%            for i=1:num
%                extend = sum(clustered_map(:)==i);
%                extend_map = extend_map + (clustered_map == i).*extend; % not a 'true' sum since there is no overlap
%            end
           tfce(:,:,index) = (extend_map.^E).*h^H.*dh;
           index = index +1;
       end
       
       % compute final score
       tfce_score = sum(tfce,3);
       
   else
       
       pos_data = (data > 0).*data;
       neg_data = abs((data < 0).*data);
       clear data
       
       % select a height, obtain cluster map, obtain extend map
       % then tfce score for that height
       l = length(min(pos_data(:)):increment:max(pos_data(:)));
       pos_increment = (max(pos_data(:)) - min(pos_data(:))) / l;
       pos_tfce = NaN(x,y,l); index = 1; 
       for h=min(pos_data(:)):pos_increment:max(pos_data(:))
           [clustered_map, num] = bwlabel(pos_data > h);
           extend_map = zeros(x,y); % same as cluster map but contains extend value instead
           extend_map = integrate(clustered_map,num,extend_map);
%            for i=1:num
%                extend = sum(clustered_map(:)==i);
%                extend_map = extend_map + (clustered_map == i).*extend;
%            end
           pos_tfce(:,:,index) = (extend_map.^E).*h^H.*dh;
           index = index +1;
       end

       l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
       neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
       neg_tfce = NaN(x,y,l); index = 1; 
       for h=min(neg_data(:)):neg_increment:max(neg_data(:))
           [clustered_map, num] = bwlabel(neg_data > h);
           extend_map = zeros(x,y); % same as cluster map but contains extend value instead
           extend_map = integrate(clustered_map,num,extend_map);
%            for i=1:num
%                extend = sum(clustered_map(:)==i);
%                extend_map = extend_map + (clustered_map == i).*extend;
%            end
           neg_tfce(:,:,index) = (extend_map.^E).*h^H.*dh;
           index = index +1;
       end
       
       % compute final score
       tfce_score = sum(pos_tfce,3)+sum(neg_tfce,3);       
   end

        
    case{2}
   % ------- tfce bootstrapped data under H0 --------------
   tfce_score = NaN(x,y,b);
   hwaitbar = waitbar(0,'TFCE');
   % check negative values if so do negate and add scores
   if min(data(:)) >= 0
       
       % select a height, obtain cluster map, obtain extend map
       % then tfce score for that height
       for boot=1:b
           tmp_data = squeeze(data(:,:,boot));
           % define increment size
           increment = (max(tmp_data(:)) - min(tmp_data(:))) / precision;
           index = 1; tfce = NaN(x,y,length(min(tmp_data(:)):increment:max(tmp_data(:))));
           %fprintf('estimating tfce under H0 boot %g \n',boot)
           
           for h=min(tmp_data(:)):increment:max(tmp_data(:))
               [clustered_map, num] = bwlabel(tmp_data > h);
               extend_map = zeros(x,y); % same as cluster map but contains extend value instead
               extend_map = integrate(clustered_map,num,extend_map);
%                for i=1:num
%                    extend = sum(clustered_map(:)==i);
%                    extend_map = extend_map + (clustered_map == i).*extend;
%                end
               tfce(:,:,index) = (extend_map.^E).*h^H.*dh;
               index = index +1;
           end
           tfce_score(:,:,boot) = sum(tfce,3);
           waitbar(boot / b);
       end
       
       
   else
       
       for boot=1:b
           %fprintf('estimating tfce under H0 for boot %g \n',boot)
           tmp_data = squeeze(data(:,:,boot));
           
           % define increment size
           increment = (max(tmp_data(:)) - min(tmp_data(:))) / precision;
           
           pos_data = (tmp_data > 0).*tmp_data;
           neg_data = abs((tmp_data < 0).*tmp_data);
           clear tmp_data
           
           % select a height, obtain cluster map, obtain extend map
           % then tfce score for that height
           l = length(min(pos_data(:)):increment:max(pos_data(:)));
           pos_increment = (max(pos_data(:)) - min(pos_data(:))) / l;
           pos_tfce = NaN(x,y,l); index = 1;
           for h=min(pos_data(:)):pos_increment:max(pos_data(:))
               [clustered_map, num] = bwlabel(pos_data > h);
               extend_map = zeros(x,y); % same as cluster map but contains extend value instead
               extend_map = integrate(clustered_map,num,extend_map);
%                for i=1:num
%                    extend = sum(clustered_map(:)==i);
%                    extend_map = extend_map + (clustered_map == i).*extend;
%                end
               pos_tfce(:,:,index) = (extend_map.^E).*h^H.*dh;
               index = index +1;
           end
           
           l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
           neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
           neg_tfce = NaN(x,y,l); index = 1;
           for h=min(neg_data(:)):neg_increment:max(neg_data(:))
               [clustered_map, num] = bwlabel(neg_data > h);
               extend_map = zeros(x,y); % same as cluster map but contains extend value instead
               extend_map = integrate(clustered_map,num,extend_map);
%                for i=1:num
%                    extend = sum(clustered_map(:)==i);
%                    extend_map = extend_map + (clustered_map == i).*extend;
%                end
               neg_tfce(:,:,index) = (extend_map.^E).*h^H.*dh;
               index = index +1;
           end
           
           % compute final score
           tfce_score(:,:,boot) = nansum(pos_tfce,3)+nansum(neg_tfce,3);
            waitbar(boot / b);
       end
   end
   close(hwaitbar)
   
end
end
%% faster integration a la Bruno Giordano
function extent_map = integrate(clustered_map,num,extent_map)

clustered_map=clustered_map(:);
nv=histc(clustered_map,0:num);
[~,idxall]=sort(clustered_map,'ascend');
idxall(1:nv(1))=[];
nv(1)=[];
ends=cumsum(nv);
inis=ends-nv+1;
for i=1:num
    idx=idxall(inis(i):ends(i));
    extent_map(idx)=nv(i);
end
end