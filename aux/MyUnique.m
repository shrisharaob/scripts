function b = MyUnique(a, varargin)
    
    
     [ IF_PRESERVE_ORDER] = DefaultArgs(varargin, {true});
     
     if IF_PRESERVE_ORDER
         [~, ii] = unique(a, 'first');
         b = a(sort(ii));
     else
         b = unique(a);
     end
end
    