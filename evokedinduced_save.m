% evokedinduced_save() - function which saves some information about calculs
%                        used with timefrequencyevoked(), timefrequencyinduced(),
%                        timefrequencyevokedBase(), timefrequecyinducedBase()
%                        or evokedinducedtopoplot(). These information are
%                        saved in a file chosen by the user. The
%                        corresponding figure is also edited to the same
%                        path as the previous file and with the same name.
%                        After saving, a dialog box appears, indicating
%                        where are saved the file and the edited figure and
%                        propose closing the figure. 
%
% Author: Fanny Grosselin, 2018
%
% Usage:
%   >> evokedinduced_save(savedata_matrix);
%
% Input:
%   savedata_matrix - [structure] structure containing figure handles, datafile names,
%                 markerfile names, headerfile names, the type of event, the sampling 
%                 rate, the number of epochs, the lower frequency in the map, the higher 
%                 frequency in the map, the frequency resolution in the map, if the user 
%                 have chosen baseline or not. If a baseline is chosen, the function 
%                 saves the type of the baseline correction 
%                 (Z-score, absolute, relative, dB) and if it's an automatic or a
%                 manual baseline. If a significance level is chosen, it will save this value.
%                 If it's an automatic baseline, the function saves the channel
%                 number used for choosing baseline, the marker of the baseline start, the sign of
%                 baseline latency, the baseline length. If it's a manual baseline, the function 
%                 saves the baseline start (in ms) and the baseline end (in ms). It also saves 
%                 the channel(s) that the user want to compute, the name of algorithm chosen and if 
%                 the user want to display significant line or not.
%                 WARNING : the figure handle indicates the location of these information 
%                 in the structure (examples : the information of figure number 2 are 
%                 saved in savedata_matrix(1,2); the information of figure number 10 are 
%                 saved in savedata_matrix(1,10)).
%
%                 Construct structure like this : 
%                 savedata.figure = handle figure; 
%                 savedata.datafile = datafile name;
%                 savedata.markerfile = markerfile name;
%                 savedata.headerfile = headerfile name;
%                 savedata.markerstimulus = EEG.event(1,2).type; %(it's an example)
%                 savedata.samplingrate = EEG.srate;
%                 savedata.epochnumber = EEG.trials;
%                 savedata.lowfreq = lower_frequency_in_Hz * EEGsrate;
%                 savedata.highfreq = higher_frequency_in_Hz * EEGsrate;
%                 savedata.resolfreq = frequency_resolution_in_Hz * EEGsrate;
%                 savedata.baseline = 'yes'; % 'automatic' or 'manual' or 'no';
%                 savedata.minb = '-700'; % number|'NaN'|''|[];
%                 savedata.maxb = '-400'; % number|'NaN'|''|[];
%                 savedata.channelchosenforbaseline = '' or [] or 'flow' or 'pressure';
%                 savedata.beginningbaseline = '' or [] or 'Start of inspiration' or 'Start 
%                                              of exhalation' or 'Maximum of 
%                                              flow/pressure' or 'Minimum of flow/pressure';
%                 savedata.latencysign = '' or []or 'before' or 'after';
%                 savedata.baselinelength = '' or [] or an integer as '300';
%				  savedata.correction_baseline = '' or 'z-score', or 'absolute', or 'relative', or 'dB'.
%                 savedata.signif_level = a float as 0.05;
%                 savedata.channelplot = channel name you wan to compute;
%                 savedata.algorithm = 'Evoked activities'; % Or 'Induced
%                       activities' or 'Evoked activities (corrected with
%                       baseline)' or 'Induced activities (corrected with
%                       baseline)'.
%                 savedata.significantline = 'yes'; % 'yes' or 'no'
%                 % recording structure in the data matrix
%                 savedata_matrix(1,handle_figure) = savedata;
%
% Output from calling:
%   The function calls a dialog box where user can choose the path and the 
%   name of the file will contain previous described information. The file
%   is saved and the corresponding figure is also edited to the same path
%   as the previous file and with the same name.
%
% See also:
%   timefrequencyevoked(), timefrequencyinduced(), timefrequencyevokedBase(),
%   timefrequecyinducedBase(), evokedinducedtopoplot().

function evokedinduced_save(savedata_matrix)

hi0 = getappdata(0,'hi');
hi = hi0.Number;

%% Save data
%  ---------
[filename,pathname] = uiputfile({'*.txt';'*.doc';'*.*'},'Save data as...');

if ischar(filename) && ischar(pathname) % to evoid error Matlab if you cancel the safeguard
    % Replace data
    % ------------
    fid = fopen([pathname,filename],'w');
    fclose(fid);

    % Write data
    % ----------
    fid = fopen([pathname,filename],'a');
    fprintf(fid,'Data file : %s\n', savedata_matrix(hi).datafile);
    fprintf(fid,'Marker file : %s\n', savedata_matrix(hi).markerfile);
    fprintf(fid,'Header file : %s\n', savedata_matrix(hi).headerfile);
    fprintf(fid,'\n');

    fprintf(fid,'Marker stimulus : %s\n',savedata_matrix(hi).markerstimulus);
    fprintf(fid,'Sampling rate : %d Hz\n',savedata_matrix(hi).samplingrate);
    fprintf(fid,'Epoch number : %d\n',savedata_matrix(hi).epochnumber);
    fprintf(fid,'\n');

    fprintf(fid,'------------------------------------------------------------------------\n');
    fprintf(fid,'Frequency parameters :\n');
    fprintf(fid,'	- Lower frequency : %.3f Hz\n',savedata_matrix(hi).lowfreq);
    fprintf(fid,' 	- Higher frequency: %.3f Hz\n',savedata_matrix(hi).highfreq);
    fprintf(fid,'	- Resolution frequency : %.3f Hz\n',savedata_matrix(hi).resolfreq);
    fprintf(fid,'Baseline : %s\n',savedata_matrix(hi).baseline);
    if strcmp(savedata_matrix(hi).baseline,'automatic')
        fprintf(fid,'	- Channel chosen for baseline choice : %s\n', savedata_matrix(hi).channelchosenforbaseline);%
        fprintf(fid,' 	- Beginning of baseline : %s\n',savedata_matrix(hi).beginningbaseline);
        fprintf(fid,'	- Latency sign : %s\n', savedata_matrix(hi).latencysign);
        fprintf(fid,'	- Baseline length in ms : %s\n', savedata_matrix(hi).baselinelength);
    elseif strcmp(savedata_matrix(hi).baseline,'manual')
        fprintf(fid,'	- Baseline start in ms : %s\n', savedata_matrix(hi).minb);
        fprintf(fid,' 	- Baseline end in ms : %s\n',savedata_matrix(hi).maxb);
    end;
    fprintf(fid,'	- Correction of the baseline : %s\n', savedata_matrix(hi).correction_baseline);
    if savedata_matrix(hi).signif_level==0
        fprintf(fid,'	- Show significance points with respect to the baseline : No\n');
    else
        fprintf(fid,'	- Show significance points with respect to the baseline : Yes\n');
        fprintf(fid,'	- Significance level : %s\n', savedata_matrix(hi).signif_level);
    end;
    fprintf(fid,'Plot :\n');
    fprintf(fid,'	- Channel plot : %s\n',savedata_matrix(hi).channelplot);
    fprintf(fid,' 	- Algorithm : %s\n', savedata_matrix(hi).algorithm);
    fprintf(fid,'	- Cone of influence : %s\n', savedata_matrix(hi).significantline);
    fprintf(fid,'------------------------------------------------------------------------\n');

    fclose(fid);

    %% Save figure
    %  -----------
    i = 1;
    filenamefig = '';
    while ~strcmp(filename(i),'.') & i < 100  % Remove extension 
        filenamefig = [filenamefig, filename(i)];
        i = i+1;
    end;
    filenamefig = [filenamefig, '.jpg'];
    saveas(hi,fullfile(pathname,filenamefig));

    %% Close current figure : dialog box
    %  ---------------------------------
    questclosefig = questdlg2(sprintf('Data are saved in %s as %s.\nFigure is saved in %s as %s.\n\nDo you want to close this figure now ?',pathname,filename,pathname,filenamefig),...
    'Close current figure','Cancel','Yes','No','Yes');

    switch questclosefig,
        case 'Yes',
            close gcf;
    end; % switch
end;

return;
