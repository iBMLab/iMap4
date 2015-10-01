function combMat = allcombs(varargin)
  varargin = varargin{1};
  sizeVec = cellfun('prodofsize', varargin);
  indices = fliplr(arrayfun(@(n) {1:n}, sizeVec));
  if isempty(indices)==0
  [indices{:}] = ndgrid(indices{:});
  combMat = cellfun(@(c,i) {reshape(c(i(:)), [], 1)}, ...
                    varargin, fliplr(indices));
  combMat = [combMat{:}];
  else
      combMat = [];
  end