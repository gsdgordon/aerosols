% Modify the string below so that it points to the parent directory of 
% your CmdStan installation
% TODO
%  o some basic checking
%  o some way to manage fileseparators?
function d = stan_home()

if ~ispc
    d = '/home/george/Apps/cmdstan';
else
    d = 'E:\Apps\cmdstan-2.24.1';
end
%d = 'C:\Users\brian\Downloads\cmdstan';
