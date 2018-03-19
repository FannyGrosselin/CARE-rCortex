% timefrequencyinducedBase() - function which calculate time-frequency maps of
%                             induced activities normalized with baseline and plot it.
%							  We compute time-frequency maps and the cardio-respiratory event
%							  potentials on around 100 trials to be faster.
%
% Author: Fanny Grosselin, 2018
%                       
% Usage:
%   >> com = timefrequencyinducedBase(marker_baseline, latpos, listauxbase, significant_line, unit, aux_plot,EEGnbchan, EEGtrials, EEGdata, resulparamclbk, nfft, EEGchanlocslabels, EEGchanlocs, listchannel, range_colorbar, EEGtimes, EEGsrate, EEGchaninfo, len_base, MINB, MAXB, minb, maxb, correction_baseline, signif_level)
%
% Inputs:
%   marker_baseline - [char] Write :
%               - 'Start of inspiration' if you want to choose baseline
%                 according to the start of inspiration. 
%            or - 'Start of exhalation' if you want to choose baseline
%                 according to the start of exhalation.
%            or - 'Maximum of flow/pressure' if you want to choose baseline
%                 according to the maximum of flow or pressure signal.
%            or - 'Minimum of flow/pressure' if you want to choose baseline
%                 according to the minimum of flow or pressure signal.
%            or - 'Cardiac' if you want to choose baseline on cardiac signal.
%   latpos [char] 'positive' : locate the baseline on the right of markers
%               which are described above. 'negative' : locate the baseline
%               on the left of markers which are descirbed above. 
%   listauxbase [char] 'pressure' : choose the baseline according to the
%               pressure signal. 'flow' : choose the baseline according to
%               the flow signal. 'cardiac' : choose the baseline according to
%               the cardiac signal.
%   significant_line - [char] 'yes' : display significant line
%               'no' : don't display significant line.
%   unit - [cell of char] Each cell contains unit of the signal 
%              of the corresponding channel (examples : first cell contains 
%              unit of the 1st number channel; the fifth cell contains unit 
%              of the 15th number channel ). Warning : write between {}.
%   aux_plot - [integer] 0 if you don't want to display, under the
%              time-frequency map, a plot of mean of channel epochs with baseline
%              location. Else, the channel number you want to plot the
%              mean through epochs. 
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
%   len_base [char] Write the baseline length in ms.
%   minb - [integer] Write the point corresponding to the baseline start. 
%          It must be between EEG.xmin*1000 and EEG.xmax*1000.
%   maxb - [integer] Write the point corresponding to the baseline end. 
%          It must be between EEG.xmin*1000 and EEG.xmax*1000.
%   MINB - [integer] Write the time point corresponding to the baseline
%          start. It must be between 0 and EEG.pnts.
%   MAXB - [integer] Write the time point corresponding to the baseline
%          end. It must be between 0 and EEG.pnts.
%   correction_baseline - [string] "z-score"|"absolute"|"relative"|"dB"
%			  "z-score" means that the mean value of baseline is substracted 
%			   from data and the standard deviation of baseline is divided
%			  "absolute" means that the mean value of baseline is subtracted
%			   from data
%		      "relative" means the data is a percentage of the mean value of baseline.
%			  "dB" means the data is divided by the mean value of baseline and then 
%			       transfomed in decibel (10*log10 operation).
%   signif_level - [value] Significance level to compute two-tailed
%           permutation significance probability level. This allows to
%           show significant points of the time-frequency map based on
%           this threshold. It will compute the significant mask based on baseline 
%           permutation, described in Grandchamp and Delorme, Frontiers in 
%           Pyscology (Sept 2011). Here we used 500 baseline permutations for 
%           each frequency bin. 
%           The surrogate time-frequency map with randome baseline is compared 
%           to the time-frequency map, normalized with a z-score baseline 
%           correction.
%           The obtained mask can be multiplied on the original
%           time-frequency map normalized by the chosen type baseline. 
%
% Outputs from pop-up:
%   com       - history string
%   Plot(s) containing a time-frequency map normalized with chosen baseline
%   and a scalp map showing channel location. If aux_plot is different of 0,
%   a plot of mean of a channel epochs with baseline location is displayed
%   thanks to minb and maxb.
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
%   pop_evokedinduced(), waveletTransformFourierPAD(), evokedinducedtopoplot(), evokedinduced_save(), plotaux().

function com = timefrequencyinducedBase( marker_baseline, latpos, listauxbase, significant_line, unit, aux_plot,EEGnbchan, EEGtrials, EEGdata, resulparamclbk, nfft, EEGchanlocslabels, EEGchanlocs, listchannel, range_colorbar, EEGtimes, EEGsrate, EEGchaninfo, len_base, MINB, MAXB, minb, maxb, correction_baseline, signif_level)
EEG = evalin('base','EEG');
ALLEEG = evalin('base','ALLEEG');
CURRENTSET = evalin('base','CURRENTSET');

for e = 1:length(EEG.chanlocs)
	EEG.CARE_rCortex(e).labels = EEG.chanlocs(e).labels;
end;

R = 200; % number of permutation for the significance of the TF map related to the choice of the baseline
nbEpoch = 100; % number or epochs chosen for the calculation of the TF maps and the epoch averaging --> to be faster

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
M = {};
if signif_level == 0
	progressbar('Computing for the chosen channel(s)...')
else
	progressbar('Computing for the chosen channel(s)...','Computing significant parts...')
end;
for k = 1:length(truelistchannel);   
    if signif_level == 0 % no need to see the significant parts of the time-frequency map
        A = 0;
        for m = 1:length(epochChosen)
			fprintf('Processing of the %s th trial...\n', num2str(m));
            [Y,freqvector] = waveletTransformFourierPAD(resulparamclbk(1).valparam,resulparamclbk(2).valparam,resulparamclbk(3).valparam,nfft,EEGdata(truelistchannel(k),:,epochChosen(m)));
            Y = (sqrt((Y) .* conj(Y))).^2; % faster than abs(Y).^2;
            % ---------------------------------------
            % Normalization of data with baseline
            % ---------------------------------------
            disp('Normalization of data with the chosen baseline, please wait...');
			% Calcul of baseline
			B = Y(:,MINB(m):MAXB(m)); % Construction of baseline vector to compute std
			% Normalization of data with z-score baseline
			mB = mean(B,2);
			stdB = std(B,0,2);
			if strcmp(correction_baseline,'z-score')
				Y = bsxfun(@rdivide, (Y-mB), stdB); % without unit       %Y(t,:) = (Y(t,:)-mean(B))./std(B,1); 
			elseif strcmp(correction_baseline,'absolute')
				Y = Y-mB; %unit^2/Hz      %Y(t,:) = Y(t,:)-mean(B);
			elseif strcmp(correction_baseline,'relative')
				Y = bsxfun(@rdivide, Y, mB).*100; % percent     % Y(t,:) = Y(t,:)./mean(B).*100;
			elseif strcmp(correction_baseline,'dB')
				Y = log10(bsxfun(@rdivide, Y, mB)).*10; %dB    %Y(t,:) = 10.*log10(Y(t,:)./mean(B))
			end;
		
            A = A + Y;
        end;
        A = A./length(epochChosen);
		EEG.CARE_rCortex(truelistchannel(k)).TFmask = [];
    else % we want to compute the significant parts of the time-frequency map
		progressbar([],0) % Update figure 
		% Baseline significance answers to this question:
		%   - Is a given Event Related Spectral Perturbation (ERSP) computed with a given baseline, significantly 
		%     different from those values obtained by a random baseline ?
		%
		% Interpretation:
		%   - If no regions in the T-F plane show significant differences, probably our baseline limits are not 
		%     pertinent or the ERSP activity has simply no relationship with the time-locked event   
		%
		% Method based on baseline permutation, described in Grandchamp and Delorme, Frontiers in Pyscology (Sept 2011) 
		
		Nf = length(resulparamclbk(1).valparam:resulparamclbk(3).valparam:resulparamclbk(2).valparam);
		Np = size(EEGdata,2);
		tBL = MAXB(1)-MINB(1)+1;
		% Construct dataset with random baseline
		NB = (tBL - 1) + 1;
		% correct signif_level for multiple comparisons (points in freq)
		signif_level = signif_level/Np;

		% Start by constructing the histograms for each freq
		% Initialise histograms
		f_hist = zeros(Nf,R); 
		
		Araw = zeros(length(epochChosen),Nf,Np);
		AdB = zeros(Nf,Np);
        A = zeros(Nf,Np);
        numEp = 1;
		
		for m = 1:length(epochChosen)
				%fprintf('Processing of the %s th trial...\n', num2str(m));
				[Y,freqvector] = waveletTransformFourierPAD(resulparamclbk(1).valparam,resulparamclbk(2).valparam,resulparamclbk(3).valparam,nfft,EEGdata(truelistchannel(k),:,epochChosen(m)));
				Y = (sqrt((Y) .* conj(Y))).^2; % faster than abs(Y).^2;
				YdB = Y;
				Ynorm = Y;
				
				% -------------------------------------------------------------
				% Normalization of data z-score baseline to prepare the inputs
				% for the baseline_significance.m function and normalization of
				% data with the chosen baseline.
				% -------------------------------------------------------------
				disp('Preparation of data to the computation of significant parts, please wait...');
				fprintf('Normalization of data with z-score...\n')
				% Calcul of baseline
				B = YdB(:,MINB(m):MAXB(m)); % Construction of baseline vector to compute std
				% Normalization of data with z-score baseline
				mB = mean(B,2);
				stdB = std(B,0,2);
				YdB = log10(bsxfun(@rdivide, YdB, mB)).*10; %dB    %YdB(t,:) = 10.*log10(YdB(t,:)./mean(B))
				
				% Normalization of data with the chosen baseline
				fprintf('Normalization of data with the chosen baseline...\n')
				if strcmp(correction_baseline,'z-score')
					Ynorm = bsxfun(@rdivide, (Ynorm-mB), stdB); % without unit       %Ynorm(t,:) = (Ynorm(t,:)-mean(B))./std(B,1); 
				elseif strcmp(correction_baseline,'absolute')
					Ynorm = Ynorm-mB; %unit^2/Hz      %Ynorm(t,:) = Ynorm(t,:)-mean(B);
				elseif strcmp(correction_baseline,'relative')
					Ynorm = bsxfun(@rdivide, Ynorm, mB).*100; % percent     % Ynorm(t,:) = Ynorm(t,:)./mean(B).*100;
				elseif strcmp(correction_baseline,'dB')
					Ynorm = YdB; % dB
				end;
				
				Araw(numEp,:,:) = Y;
				AdB = AdB + YdB;
				A = A + Ynorm;
                numEp = numEp + 1;
		end 
			
		% Random values to start baselines in EEG and random the indexes of trials
		BL_interv_ini = zeros(length(epochChosen),1);
		BL = zeros(length(epochChosen),size(Araw,2),NB);
		BL_mean = zeros(length(epochChosen),Nf);
		for i=1:R
			BL_interv_ini = randi([1 Np-NB],length(epochChosen),1)';
			BL = Araw(:,:,BL_interv_ini:BL_interv_ini+NB-1);
			BL_mean  = nanmean(BL,3);
			ind_tr = shuffle([1:length(epochChosen)]);
			
			tmp_permA = log10(bsxfun(@rdivide,Araw(ind_tr,:,:),BL_mean)).*10;
			permA = squeeze(nanmean(tmp_permA,1));
			f_hist(:,i) = nanmean(permA,2);
			progressbar([],i/R) % Update figure 
		end;
		
        AdB = AdB./length(epochChosen);
        A = A./length(epochChosen);
        % -----------------------------------------------------------
        
		
		% Find tails
		tails = quantile(f_hist,[signif_level/2 1-signif_level/2],2);

		% Find TF points in original TF map within the tails of surrogate distribution
		mask=bsxfun(@lt,AdB,tails(:,1))| bsxfun(@gt,AdB,tails(:,2));
	
		M{k} = mask;
		EEG.CARE_rCortex(truelistchannel(k)).TFmask = mask;
    end;

    Z{k} = A;
	EEG.CARE_rCortex(truelistchannel(k)).TFmap = A;
	
    channelname = trueEEGchanlocslabels(k);
    fprintf('timefrequencyinducedBase : time frequency decomposition on channel %s\n', channelname{:});
	if signif_level == 0
		if k == length(truelistchannel)
			progressbar(0.9)
			progressbar(k/length(truelistchannel))
		else
			progressbar(k/length(truelistchannel))
		end;
	else
		if k == length(truelistchannel)
			progressbar(0.9,[])
			progressbar(k/length(truelistchannel),[])
		else
			progressbar(k/length(truelistchannel),[])
		end;
	end;
end;
com = sprintf('timefrequencyinducedBase : time frequency decomposition\n');      

ALLEEG(CURRENTSET).CARE_rCortex = {};
ALLEEG(CURRENTSET) = EEG;
assignin('caller','EEG',EEG)
assignin('caller','ALLEEG',ALLEEG)
signif_level = signif_level*size(EEGdata,2);

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
    if aux_plot == 0
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
		if signif_level ~= 0
			Xmyimage = get(myimage,'XData'); 
			Ymyimage = get(myimage,'YData'); 
			I = (M{:}<1);
			% Make a truecolor all-gray image.
			greenI = cat(3, 0.5*ones(size(M{:})), 0.5*ones(size(M{:})), 0.5*ones(size(M{:})));
			hold on
			hh = imagesc(greenI);
			set(hh,'XData',Xmyimage);
			set(hh,'YData',Ymyimage);
			% Use our influence map image as the AlphaData for the solid gray image.
			set(hh, 'AlphaData', 0.8*I);
		end;
		axis('xy'); % inversion of direction y
        title(sprintf('Time frequency induced map (corrected with %s baseline) on channel %s',correction_baseline, char(EEGchanlocslabels(truelistchannel))),'edgecolor','red','linewidth',1);
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
        if ~isnan(Z{:})
            if isstr(range_colorbar)
                caxis([min(min(Z{:})) max(max(Z{:}))]);
            else
                caxis([range_colorbar(1) range_colorbar(2)]);
            end;
            c = colorbar;
            if strcmp(correction_baseline,'z-score')
				title(c,'z-score');
			elseif strcmp(correction_baseline,'absolute')
				ti = sprintf('%s^{2}/Hz', char(unit{truelistchannel}));
				title(c,ti);
			elseif strcmp(correction_baseline,'relative')
				title(c,'%');
			elseif strcmp(correction_baseline,'dB')
				title(c,'dB');
			end;
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
            if strcmp(correction_baseline,'z-score')
				title(c,'z-score');
			elseif strcmp(correction_baseline,'absolute')
				ti = sprintf('%s^{2}/Hz', char(unit{truelistchannel}));
				title(c,ti);
			elseif strcmp(correction_baseline,'relative')
				title(c,'%');
			elseif strcmp(correction_baseline,'dB')
				title(c,'dB');
			end;      
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
        if strcmp(marker_baseline,' ')
            savedata.baseline = 'manual';
            savedata.minb = num2str(minb(1));
            savedata.maxb = num2str(maxb(1));
        else
            savedata.baseline = 'automatic';
            savedata.minb = ' ';
            savedata.maxb = ' ';
        end;
		savedata.correction_baseline = correction_baseline;
        savedata.signif_level = signif_level;
        savedata.channelchosenforbaseline = listauxbase;
        savedata.beginningbaseline = marker_baseline;
        savedata.latencysign = latpos;
        savedata.baselinelength = len_base;
        savedata.channelplot = char(EEGchanlocslabels(truelistchannel));
        savedata.algorithm = 'Induced activities (corrected with baseline)';
        savedata.significantline = significant_line;
 
        % recording structure in the data matrix
        savedata_matrix(1,fig.Number) = savedata;
        setappdata(0,'savedata_matrix',savedata_matrix);
        
        % Save button
        % -----------
        uicontrol('style','pushbutton','string', 'Save data & figure','backgroundcolor',[0 1 0.5],'units','normalized','position', [4*scrsz(3)/5, scrsz(4)/70, scrsz(3)/6, scrsz(4)/15],'fontweight', 'bold','callback',clbk_savedata);

    else
		EEG = evalin('base','EEG');
        % MAP PLOT
        % --------
        set(0,'Units','points');
        scrsz = get(0,'ScreenSize');
        fig = figure('Position',[1 1 scrsz(3)/4 scrsz(4)]); 
        set(0,'Units','normalized')
        scrsz = get(0,'ScreenSize');
        h1 = subplot(2,1,1);
        myimage = imagesc(EEGtimes, freqvector*EEGsrate, Z{:});
		m=100;
		cm_plasma=plasma(m);
        colormap(h1,plasma);
		if signif_level ~= 0
			Xmyimage = get(myimage,'XData'); 
			Ymyimage = get(myimage,'YData'); 
			I = (M{:}<1);
			% Make a truecolor all-gray image.
			greenI = cat(3, 0.5*ones(size(M{:})), 0.5*ones(size(M{:})), 0.5*ones(size(M{:})));
			hold on
			hh = imagesc(greenI);
			set(hh,'XData',Xmyimage);
			set(hh,'YData',Ymyimage);
			% Use our influence map image as the AlphaData for the solid gray image.
			set(hh, 'AlphaData', 0.8*I);
		end;
		axis('xy'); % inversion of direction y
        title(sprintf('Time frequency induced map (corrected with %s baseline) on channel %s',correction_baseline, char(EEGchanlocslabels(truelistchannel))),'edgecolor','red','linewidth',1);
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
        if ~isnan(Z{:})
            if isstr(range_colorbar)
                caxis([min(min(Z{:})) max(max(Z{:}))]);
            else
                caxis([range_colorbar(1) range_colorbar(2)]);
            end;
            c = colorbar;
            if strcmp(correction_baseline,'z-score')
				title(c,'z-score');
			elseif strcmp(correction_baseline,'absolute')
				ti = sprintf('%s^{2}/Hz', char(unit{truelistchannel}));
				title(c,ti);
			elseif strcmp(correction_baseline,'relative')
				title(c,'%');
			elseif strcmp(correction_baseline,'dB')
				title(c,'dB');
			end;
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
            if strcmp(correction_baseline,'z-score')
				title(c,'z-score');
			elseif strcmp(correction_baseline,'absolute')
				ti = sprintf('%s^{2}/Hz', char(unit{truelistchannel}));
				title(c,ti);
			elseif strcmp(correction_baseline,'relative')
				title(c,'%');
			elseif strcmp(correction_baseline,'dB')
				title(c,'dB');
			end;       
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
        
        % AUXILIARY PLOT
        % --------------
        subplot(2,1,2);
        EEG = evalin('base','EEG');
        aux_plot = EEG.chanlocs(aux_plot).labels;
        plotaux(EEG,aux_plot,minb,maxb,nbEpoch);
        
        % PLOT SCALP
        % ----------
        if isempty(EEGchanlocs{truelistchannel}.X) && isempty(EEGchanlocs{truelistchannel}.Y) && isempty(EEGchanlocs{truelistchannel}.Z)
            set(gcf,'Color',[0.95,0.98,1]);
            errordlg2('Plotting scalp is not possible because there is not correct locations for this channel (X, Y and Z can not be empty)!','Warning');
        else
            set(0,'Units','normalized');
            scrsz = get(0,'ScreenSize');
            axes('Position',[scrsz(3)/1.2, scrsz(4)/2.3, scrsz(3)/8, scrsz(4)/8]);
            topoplot(truelistchannel,EEGchanlocs{truelistchannel},'electrodes','on','style', 'blank','emarker',{'.','r',20,1},'noplot','off','chaninfo',EEGchaninfo);
            title(char(EEGchanlocslabels(truelistchannel)),'color','red');
        end;
        
        % ALLOW SAVE DATA
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
        if strcmp(marker_baseline,' ')
            savedata.baseline = 'manual';
            savedata.minb = num2str(minb(1));
            savedata.maxb = num2str(maxb(1));
        else
            savedata.baseline = 'automatic';
            savedata.minb = ' ';
            savedata.maxb = ' ';
        end;
		savedata.correction_baseline = correction_baseline;
        savedata.signif_level = signif_level;
        savedata.channelchosenforbaseline = listauxbase;
        savedata.beginningbaseline = marker_baseline;
        savedata.latencysign = latpos;
        savedata.baselinelength = len_base;
        savedata.channelplot = char(EEGchanlocslabels(truelistchannel));
        savedata.algorithm = 'Induced activities (corrected with baseline)';
        savedata.significantline = significant_line;
 
        % recording structure in the data matrix
        savedata_matrix(1,fig.Number) = savedata;
        setappdata(0,'savedata_matrix',savedata_matrix);
        
        % SAVE BUTTON
        % -----------
        uicontrol('style','pushbutton','string', 'Save data & figure','backgroundcolor',[0 1 0.5],'units','normalized','position', [4*scrsz(3)/5, scrsz(4)/70, scrsz(3)/6, scrsz(4)/25],'fontweight', 'bold','callback',clbk_savedata);
    end;
 


% If 'All channels' is selected, or more than one channel are selected,
% plot scalp map with channel location which are selected. If user want
% to see time frequency map, he can choose channel(s) and push 'dispplay'
else
    figure;
    trueEEGchanlocs = [EEGchanlocs{truelistchannel(1:length(truelistchannel))}];
    trueunit = '';
    for i = 1:length(truelistchannel)
		if strcmp(correction_baseline,'z-score')
			trueunit = [trueunit; 'z-score'];
		elseif strcmp(correction_baseline,'absolute')
			ti = sprintf('%s^{2}/Hz', char(unit{1,i}{1,1}));
			trueunit = [trueunit; char(ti)];
		elseif strcmp(correction_baseline,'relative')
			trueunit = [trueunit; '%'];
		elseif strcmp(correction_baseline,'dB')
			trueunit = [trueunit; 'dB'];
		end;
    end;
    tit = sprintf('Time frequency induced map (corrected with %s baseline) on channel',correction_baseline);
    if strcmp(marker_baseline,' ')
        baseline = 'manual';
    else
        baseline = 'automatic';
    end;
    evokedinducedtopoplot([],trueEEGchanlocs, trueEEGchanlocslabels, range_colorbar, EEGtimes, freqvector, EEGsrate, resulparamclbk, Z, M, trueunit,tit,aux_plot,significant_line, baseline, marker_baseline, latpos, listauxbase, len_base, minb, maxb,correction_baseline,signif_level,'Induced activities (corrected with baseline)','chaninfo', EEGchaninfo);    
end;
end