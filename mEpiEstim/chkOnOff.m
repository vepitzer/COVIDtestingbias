function ret = chkOnOff(inp)
%CHKONOF Check for on/off input
%
%Input:
%   inp -- 'on'/'off', or 1/0 or true/false
%
%Output;
%   ret -- true if inp is 'on', 1 or true, otherwise false
%
    try
        validateattributes(inp, {'char'},{'nonempty'});
        switch lower(inp)
            case 'on'
                ret = true;
            case 'off'
                ret = false;
            otherwise
                error('Invalid option value')
        end
    catch
        try
            ret = inp;
            validateattributes(ret, {'numeric'},...
                {'>=',0,'<=',1,'integer','scalar'});
            if ret == 1
                ret = true;
            else
                ret = false;
            end
        catch
            ret = inp;
            validateattributes(ret, {'logical'},{'scalar'});
        end
    end
end
