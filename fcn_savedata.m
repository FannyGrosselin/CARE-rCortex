% fcn_savedata() - function which is called by evokedinducedtopoplot to
%                  call evokedinduced_save() in order to save figure and
%                  some information about calculs (see evokedinduced_save()). 
%
% Author: Fanny Grosselin, 2018
%
% Usage:
%   >> fcn_savedata();
%      WARNING : to call this function, a structure called savedata_matrix
%      (described in evokedinduced_save()) must exist in root (0).
%
% See also:
%   timefrequencyevoked(), timefrequencyinduced(), timefrequencyevokedBase(),
%   timefrequecyinducedBase(), evokedinducedtopoplot(), evokedinduced_save().

function fcn_savedata()

savedata_matrix = getappdata(0,'savedata_matrix'); % recovery of data matrix
hi = gcbf; % recovery of figure handle
setappdata(0,'hi',hi);
evokedinduced_save(savedata_matrix); % Call of function which save data in a file
clear hi;

return;