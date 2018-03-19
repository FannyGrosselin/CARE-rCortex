% eegplugin_CARE_rCortex() - CARE-rCortex plug-in for EEGLAB menu. 
%                      		 CARE-rCortex is the plug-in to identify physiological event 
%                      		 and to analyze evoked or induced activities
%                      		 Matlab Toolbox of Fanny GROSSELIN, 2018
%
% Usage:
%   >> eegplugin_CARE_rCortex(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.

    
    
function vers = eegplugin_CARE_rCortex(fig, trystrs, catchstrs)
    
    vers = 'CARE_rCortex';
    if nargin < 3
        error('eegplugin_CARE_rCortex requires 3 arguments');
    end;
    
    % find tools menu
    % ---------------
    menu = findobj(fig, 'tag', 'tools'); 
  
    
    % menu callback commands
    % ----------------------
    event_detection = [trystrs.no_check 'LASTCOM = pop_event_detection(EEG);' catchstrs.new_and_hist]; % a changer
    evokedinduced = [trystrs.no_check 'LASTCOM = pop_evokedinduced(EEG);' catchstrs.new_and_hist];

    
    % create menus
    % ------------
    submenu = uimenu( menu, 'Label', 'CARE-rCortex analysis', 'separator', 'on');
    uimenu( submenu, 'Label', 'Detection of physiological events'  , 'CallBack', event_detection);
    uimenu( submenu, 'Label', 'Wavelet analysis : evoked or induced activities'   , 'CallBack', evokedinduced);
   