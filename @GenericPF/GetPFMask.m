function mask = GetPFMask(rm, varargin)
    % returns a mask for pf 

   [pos, nSTD] = DefaultArgs(varargin, {[], 3});
   [nRows, nClmns] = size(rm);
   thresh = 3 * std(rm(:));
   mask = false(nRows, nClmns);
   if rm(pos) < thresh 
          return;
   end
   mask = rm > thresh;
   % right
   for x = 1 : nRows - pos(1)
       i

   

   