function varargout=aa_feval(funcname,varargin)
    try
        nout = max(nargout(funcname),1);
    catch myerr
        if (strcmp(myerr.identifier,'MATLAB:narginout:notValidMfile'))
           fprintf('%s doesn''t appear to be a valid m file?',funcname);
        else
            throw(myerr);
        end;
    end;
    % If there is a more elegant way to do this, I don't know it
    % I'd like to think there is!
    switch (nout)
        case 1
            a=feval(funcname,varargin{:});
            varargout={a};
        case 2
            [a b]=feval(funcname,varargin{:});
            varargout={a b};
        case 3
            [a b c]=feval(funcname,varargin{:});
            varargout={a b c};
        case 4
            [a b c d]=feval(funcname,varargin{:});
            varargout={a b c d};
        case 5
            [a b c d e]=feval(funcname,varargin{:});
            varargout={a b c d e};
        otherwise
            aas_log(aap,true,sprintf('Error calling aa_eval with %d arguments, function %s',nout, funcname));
    end;
end