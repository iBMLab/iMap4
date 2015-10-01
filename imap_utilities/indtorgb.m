% Junpeng Lao (2012)
% University of Fribourg
% junpeng.lao@unifr.ch

function RGB_M=indtorgb(datasetmap,cmin,cmax,datasetcmat_temp)

mapmin=min(datasetmap(:));
mapmax=max(datasetmap(:));

if isempty(cmin)==1
    cmin=mapmin;
end
if isempty(cmax)==1
    cmax=mapmax;
end

maprange=linspace(cmin,cmax,length(datasetcmat_temp));
datamaprange=maprange(maprange>=mapmin&maprange<=mapmax);
datasetcmat=datasetcmat_temp(maprange>=mapmin&maprange<=mapmax,:);

r=zeros(size(datasetmap));
g=zeros(size(datasetmap));
b=zeros(size(datasetmap));

for i=2:length(datamaprange)
    r(datasetmap>=datamaprange(i-1)&datasetmap<datamaprange(i))=datasetcmat(i-1,1);
    g(datasetmap>=datamaprange(i-1)&datasetmap<datamaprange(i))=datasetcmat(i-1,2);
    b(datasetmap>=datamaprange(i-1)&datasetmap<datamaprange(i))=datasetcmat(i-1,3);
end
r(datasetmap>datamaprange(end))=datasetcmat(end,1);
g(datasetmap>datamaprange(end))=datasetcmat(end,2);
b(datasetmap>datamaprange(end))=datasetcmat(end,3);
r(datasetmap<datamaprange(1))=datasetcmat(1,1);
g(datasetmap<datamaprange(1))=datasetcmat(1,2);
b(datasetmap<datamaprange(1))=datasetcmat(1,3);
RGB_M=zeros([size(datasetmap),3]);
RGB_M(:,:,1)=r;
RGB_M(:,:,2)=g;
RGB_M(:,:,3)=b;