function header(args)
    % A function that makes nice sub headers in the console
    % Works OK with numeric substitution, but needs a length checking
    % function.
    
    % Examples:
%         header('Header examples')
%         header({2,'Some options:'})
%         header({3,'%u %u,...',[1,2]})
    
    % Improvements; 
%             switch to varargin
%             Tests
    vals = [];
    if ischar(args)%If no header level is passed, assumed to be top-level
        msg = args;
        lvl = 0;
    elseif iscell(args)
        if ischar(args{1})
            lvl = 0;
            msg = args{1};
            p=1;
        else
            lvl = args{1};
            msg = args{2};
            p=2;
        end
        if numel(args)>p
            vals = args{p+1:end};
        end
    end

    msg_out = sprintf(msg,vals);
    blank = '                                                            ';
    marker = blank(1:3*lvl);
    fprintf([marker,' - ',msg_out,'\n'])

end