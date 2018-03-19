% pop_event_detection() - function which displays a graphical user interface
%                         allowing to choose to detect cardiac events or
%                         respiratory events on a corresponding 
%						  physiological signal.
%
% Author: Fanny Grosselin, 2018
%
% Usage:
%   >> com = pop_event_detection( EEG );
%
% Input:
%   EEG - EEGLAB EEG structure


function com = pop_event_detection(EEG)

%clear range_colorbar histo_graph event_assoc_minb event_assoc_maxb significant_line listauxbase mark marker_baseline latency listaux unit Tmin resulparamclbk evo indu lfreqnorm hfreqnorm resolnorm listchannel EEGnbchan EEGtrials EEGdata i nfft EEGtimes EEGchanlocslabels EEGchanlocs EEGsrate EEGchaninfo minb maxb MINB MAXB len_base BASELINE;

% display help if not enough arguments
% ------------------------------------
if nargin < 1
    help pop_event_detection;
    return;
    
else
    com= '';
    
    cardiac_clbk = [...
		'    wh = warndlg(''The algorithm will detect R-waves.'',''Information box'');'...
        '    uiwait(wh);'...
        '    indice_chan = 0;'...
        '    for a = 1:EEG.nbchan'...
            '    if strcmp(''CARDIAC'',EEG.chanlocs(a).labels) || strcmp(''Cardiac'',EEG.chanlocs(a).labels) || strcmp(''cardiac'',EEG.chanlocs(a).labels) ||' ...
                    'strcmp(''COEUR'',EEG.chanlocs(a).labels) || strcmp(''Coeur'',EEG.chanlocs(a).labels) || strcmp(''coeur'',EEG.chanlocs(a).labels) ||' ...
                    'strcmp(''ECG'',EEG.chanlocs(a).labels) || strcmp(''Ecg'',EEG.chanlocs(a).labels) || strcmp(''ecg'',EEG.chanlocs(a).labels) ||' ...
                    'strcmp(''CARDIAQUE'',EEG.chanlocs(a).labels) || strcmp(''Cardiaque'',EEG.chanlocs(a).labels) || strcmp(''cardiaque'',EEG.chanlocs(a).labels)'...
                    '    indice_chan = a;'...
            '    end;'...
        '    end;'...
        '    if indice_chan == 0'...
			'    str = {};'...
			'    for index = 1:EEG.nbchan'...
				'    str{end+1} = [ EEG.chanlocs(index).labels ];'...
			'    end;'...
			'    [s,v] = listdlg(''PromptString'',''No ECG detected. Choose the ECG signal:'',''SelectionMode'',''single'',''ListString'',str,''ListSize'',[300 300]);'...
			'    answer = {num2str(s)};'...
        '    else'...
			'    v = 1;'...
            '    answer = {num2str(indice_chan)};'...
        '    end;'...
        '    if isempty(answer) | v==0'...
            '    close(gcbf);'...
            '    return;'...
        '    end;'...
        '    indice_chan = str2num(cell2mat(answer));'...
        '    ecg_signal = double(EEG.data(indice_chan,:));'...
        '    [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(ecg_signal,EEG.srate,1);'...
        '    clear a answer delay ecg_signal qrs_amp_raw s wh;'...
        ... % Plot to control cardiac events
        ... % ------------------------------
        '    CARE_eegplot (EEG.data(indice_chan,:),''srate'',EEG.srate,''title'',''Remove, correct, add or accept cardiac marks -- CARE_eegplot()'',''butlabel'',''ACCEPT'',''command'',''physio_mark'',''physio_mark'',[qrs_i_raw],''physio_name'',''cardiac_event'');'];

    respiratory_clbk = [...
		'	 choice_respiratory = questdlg2(''On what type of channel do you want to detect respiratory events?'',''Choice respiratory channel'',''Pressure'',''Flow'',''Pressure'');'...
        '    indice_chan= 0;'...
		'	 switch choice_respiratory,'...
			'	case ''Pressure'','...
				'    for a = 1:EEG.nbchan'...
					'	if strcmp(''Pression'',EEG.chanlocs(a).labels) || strcmp(''pression'',EEG.chanlocs(a).labels) || strcmp(''PRESSION'',EEG.chanlocs(a).labels) ||'...
                            'strcmp(''Pressure'',EEG.chanlocs(a).labels) || strcmp(''pressure'',EEG.chanlocs(a).labels) || strcmp(''PRESSURE'',EEG.chanlocs(a).labels)'...
							'    indice_chan = a;'...
					'   end;'...
				'	end;'...
				'    if indice_chan == 0'...
					'    str = {};'...
					'    for index = 1:EEG.nbchan'...
						'    str{end+1} = [ EEG.chanlocs(index).labels ];'...
					'    end;'...
					'    [s,v] = listdlg(''PromptString'',''No pressure respiratory signal detected. Choose the pressure signal:'',''SelectionMode'',''single'',''ListString'',str,''ListSize'',[400 300]);'...
					'    answer = {num2str(s)};'...
				'    else'...
					'    v = 1;'...
					'    answer = {num2str(indice_chan)};'...
				'    end;'...
				'    if isempty(answer) | v==0'...
					'    close(gcbf);'...
					'    return;'...
				'    end;'...
				'   if str2num(cell2mat(answer)) > EEG.nbchan'...
					'   errordlg2(''Enter a correct channel number for the pressure respiratory signal.'');'...
					'   return;'...
				'   end;'...
			'	case ''Flow'','...
				'   for a = 1:EEG.nbchan'...
					'    if strcmp(''Debit'',EEG.chanlocs(a).labels) || strcmp(''debit'',EEG.chanlocs(a).labels) || strcmp(''D?bit'',EEG.chanlocs(a).labels) ||' ...
                            'strcmp(''d?bit'',EEG.chanlocs(a).labels) || strcmp(''DEBIT'',EEG.chanlocs(a).labels) || strcmp(''Flow'',EEG.chanlocs(a).labels) ||' ...
                            'strcmp(''flow'',EEG.chanlocs(a).labels) || strcmp(''FLOW'',EEG.chanlocs(a).labels) || strcmp(''D�bit'',EEG.chanlocs(a).labels) ||'...
                            'strcmp(''d�bit'',EEG.chanlocs(a).labels) || strcmp(''débit'',EEG.chanlocs(a).labels) || strcmp(''Débit'',EEG.chanlocs(a).labels)'...
							'    indice_chan = a;'...
					'   end;'...
				'	end;'...
				'    if indice_chan == 0'...
					'    str = {};'...
					'    for index = 1:EEG.nbchan'...
						'    str{end+1} = [ EEG.chanlocs(index).labels ];'...
					'    end;'...
					'    [s,v] = listdlg(''PromptString'',''No flow respiratory signal detected. Choose the flow signal:'',''SelectionMode'',''single'',''ListString'',str,''ListSize'',[400 300]);'...
					'    answer = {num2str(s)};'...
				'    else'...
					'    v = 1;'...
					'    answer = {num2str(indice_chan)};'...
				'    end;'...
				'    if isempty(answer) | v==0'...
					'    close(gcbf);'...
					'    return;'...
				'    end;'...
				'   if str2num(cell2mat(answer)) > EEG.nbchan'...
					'   errordlg2(''Enter a correct channel number for the flow respiratory signal.'');'...
					'   return;'...
				'   end;'...
        '    end;'...
        '    indice_chan = str2num(cell2mat(answer));'...
        '    respiratory_signal = double(EEG.data(indice_chan,:));'...
        '    AllChoices = {''Minimum'',''Maximum'', ''Start of inspiration'', ''Start of expiration''};'...
        '    [ChoiceNum,v]  = listdlg(''PromptString'',''What you want to detect?'','...
                ' ''SelectionMode'',''single'','...
                ' ''ListString'',AllChoices);'...
        '    if ChoiceNum == 1'...
            '    mark = ''peakmin'';'...
        '    elseif ChoiceNum == 2'...
            '    mark = ''peakmax'';'...
        '    elseif ChoiceNum == 3'...
            '    mark = ''startinspi'';'...
        '    elseif ChoiceNum == 4'...
            '    mark = ''startexhal'';'...
        '    end;'...
        '    if ~exist(''mark'')'...
            '    close(gcbf);'...
            '    return;'...
        '    end;'...
        '    [lc,vr] = find_resp_marks(respiratory_signal,EEG.srate,mark);'...
        '    clear a answer respiratory_signal AllChoices ChoiceNum s v mark vr choice_respiratory;'...
        ... % Plot to control respiratory events
        ... % ----------------------------------
        '    CARE_eegplot (EEG.data(indice_chan,:),''srate'',EEG.srate,''title'',''Remove, correct, add or accept muscular marks -- CARE_eegplot()'',''butlabel'',''ACCEPT'',''command'',''physio_mark'',''physio_mark'',[lc],''physio_name'',''resp_event'');'];

    
    
    
    %% Main GUI
    %  --------
	spg = figure;
    supergui('fig',spg,'geomhoriz', {[40] [1] [40] [1] [40] [1] [40 2 40]}, 'geomvert', [1 1 1 1 1 1 1],'uilist', { ...
        ...
        { 'Style', 'text', 'string', 'Choose a type of physiological events : ', 'fontweight', 'bold', 'fontsize',12 }...
        {}...
        { 'style', 'pushbutton' , 'string', 'Cardiac events','fontweight', 'bold','callback',cardiac_clbk}...
        {}...
        { 'style', 'pushbutton' , 'string', 'Respiratory events','fontweight', 'bold','callback',respiratory_clbk}...
        {}...
        { 'style', 'pushbutton' , 'string', 'Help' 'callback' 'pophelp(''pop_event_detection'');' }...
        {}...
        { 'style', 'pushbutton' , 'string', 'Cancel' 'callback' 'close(gcbf); return;'}...
        },'title', 'Detection of physiological events -- pop_event_detection()');
		
	set(spg,'WindowStyle','modal');
end;


return;
