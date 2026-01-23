function s = validateEnum(val, allowed, caller, argName)
%TFL.VALIDATEENUM  Lowercase + validatestring with clean error messages.

    if nargin < 3 || isempty(caller),  caller  = mfilename; end
    if nargin < 4 || isempty(argName), argName = 'option';  end

    if ~(ischar(val) || isstring(val))
        error('%s:BadOptionValue', caller, ...
              '%s must be a string. Allowed: %s', argName, strjoin(allowed, ', '));
    end

    s = lower(char(val));
    s = validatestring(s, allowed, caller, argName);
    s = lower(char(s));
end