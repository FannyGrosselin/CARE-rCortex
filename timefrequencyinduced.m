% timefrequencyinduced() - function which calculate time-frequency maps of
%                          induced activities and plot it.
%					       We compute time-frequency maps and the cardio-respiratory event
%						   potentials on around 100 trials to be faster.
%
% Author: Fanny Grosselin, 2018
%                       
% Usage:
%   >> com = timefrequencyinduced( significant_line, unit,EEGnbchan, EEGtrials, EEGdata, resulparamclbk, nfft, EEGchanlocslabels, EEGchanlocs, listchannel, range_colorbar, EEGtimes, EEGsrate, EEGchaninfo);
%
% Inputs:
%   significant_line - [char] 'yes' : display significant line
%               'no' : don't display significant line.
%   unit - [cell of char] Each cell contains unit of the signal 
%              of the corresponding channel (examples : first cell contains 
%              unit of the 1st number channel; the fifth cell contains unit 
%              of the 15th number channel ). Warning : write between {}.
%   EEGnbchan - [integer] number of channels in the dataset - see EEG.nbchan.
%   EEGtrials - [integer] number of epochs in the dataset - see EEG.trials.
%   EEGdata - [array of integers] array containing points values of each
%              epochs, of each channels, at each time - see EEG.data.
%   resulparamclbk - [structure of integers] :
%              resulparamclbk(1).valparam = lower frequency
%              resulparamclbk(2).valparam = higher frequency
%              resulparamclbk(3).valparam = frequency resolution
%              WARNING : divide your integers (in Hz) by sampling frequency
%              (see EEG.srate).
%   nfft - power of two greater than the number of points (see EEG.pnts).
%              It's an input of waveletTransformFourierPAD() function.
%   EEGchanlocslabels - [cells of char] cells containing label of each
%              channel (examples : first cell contains the name of 1st
%              channel; the fifth cell contains the name of 5th channel).
%              It's :
%              for i = 1:EEG.nbchan
%                   EEGchanlocslabels{i} = EEG.chanlocs(i).labels;
%              end;
%   EEGchanlocs - [cells of structure] cells containing information about
%              channel location (examples : first cell contains information
%              about channel number 1; fifth cell contains information
%              about channel ,number 5).
%              It's :
%              for i = 1:EEG.nbchan
%                   EEGchanlocs{i} = EEG.chanlocs(i);
%              end;
%   listchannel - [integer] 1 if you want to display the calcul on all the
%              channels. Else, write the channel number you want to
%              analyse with a shift of one (examples : if you want to 
%              analyse the channel number 1, write 2; if you want to
%              analyse the channel number 10, write 11).
%   range_colorbar - ['automatic' or array of integer] Write 'automatic' if
%              you want colorbar spreads between minimum value of time-frequency
%              map and maximum value of time-frequency map. Otherwise,
%              write an array of integer like this : [colorbar_min colorbar_max].
%   EEGtimes - [array of double] It's the same as EEG.times.
%   EEGsrate - [integer] It's the same as EEG.srate.
%   EEGchaninfo - [structure] It's the same as EEG.chaninfo. 
%
% Outputs from pop-up:
%   com       - history string
%   Plot(s) containing a time-frequency map and a scalp map showing channel  
%   location.
%   If significant_line = "yes", a line is displayed on time-frequency map 
%   separating significant data and non significant data.
%   If listchannel = 1, evokedinducedtopoplot() is called.
%   There is also a button to save figure and some data (see
%   evokedinduced_save()).
%   The function also saves a structure called savedata_matrix which
%   contains the figure handle, the datafile name, the markerfile name, the
%   headerfile name, the type of event, the sampling rate, the number of epochs, 
%   the lower frequency in the map, the higher frequency in the map,
%   the frequency resolution in the map, if the user have chosen baseline or not.
%	If a baseline is chosen, the function saves the type of the baseline correction 
%	(Z-score, absolute, relative, dB) and if it's an automatic or a
%   manual baseline. If a significance level is chosen, it will save this value.
%   If it's an automatic baseline, the function saves the channel
%   number used for choosing baseline, the marker of the baseline start, the sign of
%   baseline latency, the baseline length. If it's a manual baseline, the function 
%   saves the baseline start (in ms) and the baseline end (in ms). It also saves 
%   the channel(s) that the user want to compute, the name of algorithm chosen and if 
%   the user want to display significant line or not. WARNING : the figure handle
%   indicates the location of these information in the structure (examples
%   : the information of figure number 2 are saved in savedata_matrix(1,2);
%   the information of figure number 10 are saved in savedata_matrix(1,10)).
%
% See also:
%   pop_evokedinduced(), waveletTransformFourierPAD(), evokedinducedtopoplot(), evokedinduced_save().

function com = timefrequencyinduced( significant_line,unit, EEGnbchan, EEGtrials, EEGdata, resulparamclbk, nfft, EEGchanlocslabels, EEGchanlocs, listchannel, range_colorbar, EEGtimes, EEGsrate, EEGchaninfo)

EEG = evalin('base','EEG');
ALLEEG = evalin('base','ALLEEG');
CURRENTSET = evalin('base','CURRENTSET');

for e = 1:length(EEG.chanlocs)
	EEG.CARE_rCortex(e).labels = EEG.chanlocs(e).labels;
end;

nbEpoch = 100; % number or epochs chosen for the calculation of the TF maps --> to be faster

if EEGtrials<100
	epochChosen = [1:EEGtrials];
	fprintf('All the %s trials are chosen for the computation of the time-frequency map(s).\n',num2str(EEGtrials));
else
	epochChosen = [1:floor(EEGtrials/nbEpoch):EEGtrials]; % chose about 100 epochs by a regular way in order to don't have time bias
	fprintf('%s over %s trials are chosen for the computation of the time-frequency map(s).\n',num2str(length(epochChosen)),num2str(EEGtrials));
end;

%% Calcul of time frequency transform
%  ----------------------------------

if listchannel(1) == 2 % if 'All channels' is selected
    truelistchannel = [1:EEGnbchan];
    trueEEGchanlocslabels = [EEGchanlocslabels(1:EEGnbchan)];
else
    % if user selected one or several channel(s)
    truelistchannel = [listchannel(1:length(listchannel))-2]; % because there is a shift due to 'All channels' and ' '
    trueEEGchanlocslabels = [EEGchanlocslabels(truelistchannel(1:length(listchannel)))];
end;

disp('Please wait...');
Z = {};
progressbar('Computing for the chosen channel(s)...')
for k = 1:length(truelistchannel);
    A = 0;
    for m = 1:length(epochChosen)
        [Y,freqvector] = waveletTransformFourierPAD(resulparamclbk(1).valparam,resulparamclbk(2).valparam,resulparamclbk(3).valparam,nfft,EEGdata(truelistchannel(k),:,epochChosen(m)));
        Y = (sqrt((Y) .* conj(Y))).^2; % faster than abs(Y).^2;
        A = A + Y;
    end;
    A = A./length(epochChosen);
    Z{k} = A;
    channelname = trueEEGchanlocslabels(k);
    fprintf('timefrequencyinduced : time frequency decomposition on channel %s\n', channelname{:});  
	
	if k == length(truelistchannel)
		progressbar(0.9)
		progressbar(k/length(truelistchannel))
	else
		progressbar(k/length(truelistchannel))
	end;
	
	EEG.CARE_rCortex(truelistchannel(k)).TFmask = [];
	EEG.CARE_rCortex(truelistchannel(k)).TFmap = A;
end;
com = sprintf('timefrequencyinduced : time frequency decomposition\n');  

ALLEEG(CURRENTSET).CARE_rCortex = {};
ALLEEG(CURRENTSET) = EEG;
assignin('caller','EEG',EEG)
assignin('caller','ALLEEG',ALLEEG)
    
        
%% Callback to save data
%  ---------------------
clbk_savedata = [...
    'savedata_matrix = getappdata(0,''savedata_matrix'');'... % recovery of data matrix
    'hi = gcbf;'... % recovery of figure handle
    'setappdata(0,''hi'',hi);'...
    'evokedinduced_save(savedata_matrix);'... % Call of function which save data in a file
    'clear hi;'...
    ];

%% Plots of results
%  ----------------
% If only one channel selected, plot image the time frequency map with the
% corresponding scalp
if (length(truelistchannel) == 1)
    set(0,'Units','points');
    scrsz = get(0,'ScreenSize');
    fig = figure('Position',[1 1 scrsz(3)/4 1.2*scrsz(4)/2]); 
    set(0,'Units','normalized')
    scrsz = get(0,'ScreenSize');
    h0 = axes('Position',[scrsz(3)/5, scrsz(4)/8, scrsz(3)/1.3, scrsz(4)/1.4]);
    myimage = imagesc(EEGtimes, freqvector*EEGsrate, Z{:});
	m=100;
	cm_plasma=plasma(m);
	colormap(h0,plasma);
    axis('xy'); % inversion of direction y
    title(sprintf('Time frequency induced map on channel %s', char(EEGchanlocslabels(truelistchannel))),'edgecolor','red','linewidth',1);
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    if ~isnan(Z{:})
        if isstr(range_colorbar)
            caxis([min(min(Z{:})) max(max(Z{:}))]);
        else
            caxis([range_colorbar(1) range_colorbar(2)]);
        end;
        c = colorbar;
        ti = sprintf('%s^{2}/Hz', char(unit{truelistchannel}));
        title(c, ti);
    elseif isnan(Z{:})
        errordlg2(sprintf('WARNING : The %s signal is null or empty ! The colorbar can not be printed.',char(EEGchanlocslabels(truelistchannel))),'Warning');
    else
        Y = Z;
        for ji = 1:size(Z,1)
            for ki = 1:size(Z,2)
                if isnan(Z{ji,ki})
                    Y(ji,ki) = [];
                end;
            end;
        end;
        if isstr(range_colorbar)
            caxis([min(min(Y{:})) max(max(Y{:}))]);
        else
            caxis([range_colorbar(1) range_colorbar(2)]);
        end;
        c = colorbar;
        ti = sprintf('%s^{2}/Hz', char(unit{truelistchannel}));
        title(c, ti);
    end;

    % Plot signficant line
    % --------------------
    if strcmp(significant_line,'yes')
        ko = 5;  % for Morlet Wavelet
        fourier_factor = (4*pi)/(ko + sqrt(2+ko^2));
        coi = fourier_factor/sqrt(2);
        n1 = size(Z{:},2);
        coi = coi*1*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];
        upCurve = EEGsrate./coi;
        lowCurve = -0.5.*ones(size(coi));
        % Overlay the cone of influence
        hold on, [ph,msg]=jbfill(EEGtimes,lowCurve,upCurve,[0.2 0.2 0.2],[0 0 0],1,0.7);
        hold off;
    end;

    % Plot line on image
    % ------------------
    hold on;
    % define points (in matrix coordinates)
    p1 = [0,0];
    p2 = [min(ylim),max(ylim)];
    % plot the line
    plot([p1(1),p1(2)],[p2(1),p2(2)],'w--');
    hold off;

    % Plot scalp
    % ----------
    if isempty(EEGchanlocs{truelistchannel}.X) && isempty(EEGchanlocs{truelistchannel}.Y) && isempty(EEGchanlocs{truelistchannel}.Z)
        set(gcf,'Color',[0.95,0.98,1]);
        errordlg2('Plotting scalp is not possible because there is not correct locations for this channel (X, Y and Z can not be empty)!','Warning');
    else
        set(0,'Units','normalized')
        scrsz = get(0,'ScreenSize');
        axes('Position',[0, 3.8*(scrsz(4)/5), scrsz(3)/6, scrsz(4)/6]);
        topoplot(truelistchannel,EEGchanlocs{truelistchannel},'electrodes','on','style', 'blank','emarker',{'.','r',20,1},'noplot','off','chaninfo',EEGchaninfo);
        title(char(EEGchanlocslabels(truelistchannel)),'color','red');
    end;

    % Allow save data
    % ---------------
    EEG = evalin('base','EEG');

    % recovery of data matrix
    try
        savedata_matrix = getappdata(0,'savedata_matrix');
    end;
    if isempty(savedata_matrix)
        savedata_matrix = struct('figure',[],'datafile',[],'markerfile',[],'headerfile',[],'markerstimulus',[],'samplingrate',[],'epochnumber',[],'lowfreq',[],'highfreq',[],'resolfreq',[],'baseline',[],'minb',[],'maxb',[],'correction_baseline',[],'signif_level',[],'channelchosenforbaseline',[],'beginningbaseline',[],'latencysign',[],'baselinelength',[],'channelplot',[],'algorithm',[],'significantline',[]);
    end;

    % recovery of some data
    p = EEG.filepath;
    [linecomments colcomments] = size(EEG.comments);
    for s = 1:linecomments
		if strcmp('Original file: ',EEG.comments(s,1:15));
			datafile = deblank(EEG.comments(s,16:end));
			markerfile = [datafile(1:end-3),'vmrk'];
			headerfile = [datafile(1:end-3),'vhdr'];
		elseif strcmp('Parent dataset: ',EEG.comments(s,1:16));
			datafile = deblank(EEG.comments(s,17:end));
			markerfile = [];
			headerfile = [];
		end;
	end;

    % recording of data in a structure 
    savedata.figure = fig.Number;
    savedata.datafile = char(datafile);
    savedata.markerfile = char(markerfile);
    savedata.headerfile = char(headerfile);
    savedata.markerstimulus = EEG.event(1,2).type;
    savedata.samplingrate = EEGsrate;
    savedata.epochnumber = EEGtrials;
    savedata.lowfreq = resulparamclbk(1).valparam * EEGsrate;
    savedata.highfreq = resulparamclbk(2).valparam * EEGsrate;
    savedata.resolfreq = resulparamclbk(3).valparam * EEGsrate;
    savedata.baseline = 'no';
    savedata.minb = [];
    savedata.maxb = [];
	savedata.correction_baseline = '';
    savedata.signif_level = 0;
    savedata.channelchosenforbaseline = [];
    savedata.beginningbaseline = [];
    savedata.latencysign = [];
    savedata.baselinelength = [];
    savedata.channelplot = char(EEGchanlocslabels(truelistchannel));
    savedata.algorithm = 'Induced activities';
    savedata.significantline = significant_line;

    % recording structure in the data matrix
    savedata_matrix(1,fig.Number) = savedata;
    setappdata(0,'savedata_matrix',savedata_matrix);

    % Save button
    % -----------
    uicontrol('style','pushbutton','string', 'Save data & figure','backgroundcolor',[0 1 0.5],'units','normalized','position', [4*scrsz(3)/5, scrsz(4)/70, scrsz(3)/6, scrsz(4)/15],'fontweight', 'bold','callback',clbk_savedata);

% If 'All channels' is selected, or more than one channel are selected,
% plot scalp map with channel location which are selected. If user want
% to see time frequency map, he can choose channel(s) and push 'dispplay'
else
    figure;
    trueEEGchanlocs = [EEGchanlocs{truelistchannel(1:length(truelistchannel))}];
    trueunit = '';
    for ja = 1:length(truelistchannel)
        len_unit(ja) = length(char(unit{1,truelistchannel(ja)}));
    end;
    max_len_unit = max(len_unit);
    for i = 1:length(truelistchannel)
        if len_unit(i) ~= max_len_unit % to adapt length of unit
            for jou = 1:(max_len_unit - len_unit(i))
                unit{1,truelistchannel(i)} = [' ', unit{1,truelistchannel(i)}];
            end;
        end;
        trueunit = [trueunit; char(sprintf('%s^{2}/Hz', char(unit{1,truelistchannel(i)})))];
    end;
    tit = sprintf('Time frequency induced map on channel');
    baseline = 'no';
    aux_plot = 0;
    marker_baseline = '';
    latpos = '';    
    listauxbase = '';
    len_base = '';
	correction_baseline = '';
    signif_level = 0;
	M = {};
    evokedinducedtopoplot([],trueEEGchanlocs, trueEEGchanlocslabels, range_colorbar, EEGtimes, freqvector, EEGsrate, resulparamclbk, Z, M,trueunit, tit,aux_plot,significant_line, baseline, marker_baseline, latpos, listauxbase,len_base,'NaN','NaN',correction_baseline,signif_level,'Induced activities','chaninfo', EEGchaninfo);   
end;
end

