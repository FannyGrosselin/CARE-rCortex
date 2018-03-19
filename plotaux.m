% plotaux() - function which computes average of a channel epochs and plots it
%             with a baseline. 
%
% Author: Fanny Grosselin, 2018
%
% Usage:
%   >> plotaux(EEG,aux_plot,minb,maxb);
%
% Inputs:
%   EEG - EEGLAB EEG structure
%   aux_plot - [integer] Write the channel number you want to compute mean 
%              of epochs and plot it.
%   minb - [integer] Write the point corresponding to the baseline start. 
%          It must be between EEG.xmin*1000 and EEG.xmax*1000.
%   maxb - [integer] Write the point corresponding to the baseline end. 
%          It must be between EEG.xmin*1000 and EEG.xmax*1000.
%	nbEpoch - [integre] Write the number of epochs we will keep for the 
%			  computation and the display of event potentials.
%
% See also:
%   timefrequencyevoked(), timefrequencyinduced(), timefrequencyevokedBase(),
%   timefrequecyinducedBase(), evokedinducedtopoplot().

function plotaux(EEG,aux_plot,minb,maxb, nbEpoch)

if EEG.trials<100
	epochChosen = [1:EEG.trials];
	fprintf('All the %s trials are chosen for the computation of the time-frequency map(s).\n',num2str(EEG.trials));
else
	epochChosen = [1:floor(EEG.trials/nbEpoch):EEG.trials]; % chose about 100 epochs by a regular way in order to don't have time bias
	fprintf('%s over %s trials are chosen for the computation of the mean of epochs.\n',num2str(length(epochChosen)),num2str(EEG.trials));
end;

% Calcul of average of chosen auxiliary channel
% ---------------------------------------------
A = 0;
for i = 1:length(EEG.chanlocs)
    if strcmp(EEG.chanlocs(i).labels,aux_plot);
        a = i;
        A = sum(EEG.data(a,:,epochChosen)-mean(EEG.data(a,:,epochChosen)),3);
        A = A./length(epochChosen);
        % --------------------------------------
        % Load units
        % --------------------------------------
        p = EEG.filepath;
        [linecomments colcomments] = size(EEG.comments);
        for s = 1:linecomments
            unit = {};
            if strcmp('Original file: ',EEG.comments(s,1:15));
                o = deblank(EEG.comments(s,16:end));
                n = [o(1:end-3),'vhdr'];
                try
                    hdrunit = readbvconf(p,n);
                    if isfield(hdrunit,'channelinfos')
                        [a,b,c,u] = strread(hdrunit.channelinfos{1,a},'%s%s%s%s',1,'delimiter',',');
                    else
                        u = 'Unit';
                    end;
                    v = unicode2native(char(u));
                    if length(v) == 3 % if it is microvolts
                        unit = 'muV';
                    else
                        unit = u;
                    end;
                catch
                    unit = 'Unit';
                end;
            else
                unit = 'Unit';
            end;
        end;
    end;
end;

% Plot auxiliary channels
% -----------------------
set(0,'Units','normalized')
scrsz = get(0,'ScreenSize');
%axes('Position',[scrsz(3)/7.6, scrsz(4)/11, scrsz(3)/1.455, scrsz(4)/3]);
plot(EEG.times,A);
axis('xy'); % inversion of direction y
xlabel('Times (ms)');
ylabel(sprintf('%s %s','Intensity (',char(unit),')'));
title(sprintf('Graph of %s (average based on %s trials)',aux_plot, num2str(length(epochChosen))));
c = colorbar;
set(c,'visible','off');

hold on;

% Plot baseline on image
% ----------------------
baseline_len = (maxb(1) - minb(1));
baseline_centre = minb + baseline_len/2;
baseline_centre = round(median(baseline_centre));
ix = [round(baseline_centre - baseline_len/2) round(baseline_centre - baseline_len/2) round(baseline_centre + baseline_len/2) round(baseline_centre + baseline_len/2)];
igrec = [min(ylim) max(ylim) max(ylim) min(ylim)];
zed = [0 0 0 0];
patch(ix,igrec,zed,'y');
alpha(0.5);

% Plot line on image
% ------------------
% define points (in matrix coordinates)
p1 = [0,0];
p2 = [min(ylim),max(ylim)];
% plot the line
plot([p1(1),p1(2)],[p2(1),p2(2)],'k--');

hold off;

clear p n o hdrunit a b c u ix igrec zed;

return;