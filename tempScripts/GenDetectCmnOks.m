function funcHdl = GenDetectCmnPks(n, varargin)
% generates func hdl with variable no pf input args for detecting common peaks across sessions
    [funcName, tolerence] = DefaultArgs(varargin, {'ismemberf', 2});
    str = [];
    for i = 1 : n    
        if i == 1, str = 'A1'; 
        else
            str = [str, ',A' num2str(i)];
        end
    end
    funchdl = str2func(['@(' str ')' funcName '(' str ',' , '''row''', ',', '''tol''',  ',', num2str(tolerence), ')']);
        
end