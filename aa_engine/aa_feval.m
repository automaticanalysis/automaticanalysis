function varargout=aa_feval(funcname,varargin)
    try
        nout = max(nargout(funcname),1);
    catch myerr
        if (strcmp(myerr.identifier,'MATLAB:narginout:notValidMfile'))
           aas_log([],false,sprintf('%s doesn''t appear to be a valid m file?',funcname));
        else
            throw(myerr);
        end;
    end;

    varargout = cell(1,nout);
    [varargout{:}]=feval(funcname,varargin{:});
end