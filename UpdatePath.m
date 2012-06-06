function UpdatePath(p)
% function UpdatePath(p)
% this function add sub folders of the input path to the matlab path.
% the default value of path is current path
% input:
%       p,(optional) string that contain the absolut path of folder that its
%       sub-folder should add to the matlab path.
%       
% 
% Author: A.Rahim Kadkhodamohammadi r.k.mohammadi@gmail.com
% 07/02/2012 CBA, Uppsala University
%--------------------------------------------------------------


if nargin ==0
    path = fullfile(pwd, 'libs');
elseif nargin ==1
    path = p; 
elseif nargin > 1
    error('this function only accept on input arguemt that is library path');
end

if exist(path, 'dir') ~= 7
    error(['The input Path does not correct. path = ', path]);
end
folders = dir(path);
if numel(folders) < 3
    error('please check the path');
else
    for i=3 : numel(folders)
        if folders(i).isdir
            addpath(fullfile(path, folders(i).name));
        end
    end
end

% add libsvm library
addpath(fullfile(path,'libsvm-3.11', 'matlab'));
display('Library is updated! ');