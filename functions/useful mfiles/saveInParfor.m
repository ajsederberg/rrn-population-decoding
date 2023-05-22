function saveInParfor(varargin)
% save function that can be used in a parfor loop. Also useful in scripts
% where verifying that a variable was saved is useful. Edited on 1/9/18 to
% make '-v7.3' the default if variables are > 1 GB. 

savefile = varargin{1}; % first input argument
for i=2:nargin
    savevar.(inputname(i)) = varargin{i}; % other input arguments
end

s = whos('-var', 'savevar');
fileBytes = s.bytes;
if fileBytes > 1e8
    save(savefile,'-struct','savevar', '-v7.3')
else
    save(savefile,'-struct','savevar', '-v7')
end