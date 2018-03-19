
function [lc,vr] = find_resp_marks(pres,Fs,mark)

% Check amplitude differences due to different records in 10 sec windows
% Inputs
% pres : breathing signal (pression / flow )
% Fs : Sampling frequency
% mark : event that you want to detect. It can be 'startinspi' if you want
%       to detect the start of inspiration. It can be 'startexhal' if you want to
%       detect the start of exhalation. It can be 'peakmax' if you want to detect
%       maximum of flow/pression. It can be 'peakmin' if you want to detect
%       minimum of flow/pression.
% Outputs
% lc : It's equal to lc1 (maximum of flow/pression) or lc2 (minimum of flow/pression)
%       or lc3 (start of exhalation (sample number)) or lc4 (start of inspiration (sample
%       number)). It's an array of the sample points which are what you want to
%       detect.
% vr : It's an array of the differences between two successive points which
%      are what you want to detect.

% figure;
% plot(pres);

% Author: Xavier Navarro & Fanny Grosselin,2015

disp('Please wait, we are detecting the events...');
progressbar('Detecting respiratory events...')
					
w=10*Fs;
pres = (pres-mean(pres))/max(pres);
delta_ant =max(pres(1:w))-min(pres(1:w));
i=w; k=2;
%% Segmentation of flow/pression according to amplitude
brkpoint(1)=1;
while i<length(pres)-w
    prs = pres(i:i+w);
    delta = max(prs)-min(prs);
    ratio = abs(delta/delta_ant);
    if ratio > 3 || ratio < 0.3
        % Find segmentation point
        %disp(['Segmentation point ' num2str(k)])
        mean_E = sum(prs.^2)/Fs;
        i2=1;
        for j=1:Fs:w-Fs
            en(i2)=sum(prs(j:j+Fs).^2);
            i2=i2+1;
        end
        [~,mind] = max(en);
        %disp(['M ind ' num2str(mind)])
        [~,mnp] = min(prs(1+(mind-1)*Fs:mind*Fs));
        %disp(['Min value ' num2str(mnp)])
        brkpoint(k)=mnp+i;
        k=k+1;
        i = i+mnp;
    else
        i = i+w;
    end
    delta_ant = delta;
end
progressbar(0.1)
brkpoint(k)=length(pres);
%disp([ num2str(k) ' breakpoints'])

%% Calcul of MINPEAKDISTANCE
L = length(pres); % length of signal
f1 = 2; % cut-off frequency
Wn = 2.*f1./Fs; % normaization of this frequency
[b,a] = butter(7,Wn,'low'); % calcul of filter coefficients
signal = filter(b,a,pres); % filter
progressbar(0.2)

% FFT
y = abs(fft(signal));      %Retain Magnitude (imaginary part = phase)
y = y(1:L/2);     %Discard Half of Points
f = Fs*(0:L/2-1)/L;

% Search main frequency
[val_max,ind_max]=max(y(find(f>0.1)));
MF = f(ind_max);

if MF ~= 0
    period = 1/MF;
    MINPEAKDISTANCE = round(period*Fs/2);
else
    MINPEAKDISTANCE = 2 * Fs;
end;
%disp([ 'Main frequency = ' num2str(MF) ' Hz  ---  MINPEAKDISTANCE = ' num2str(MINPEAKDISTANCE) ' samples ']);


%% Cycle detection in segments
pk1 = []; pk3 = []; pk4 = []; pk2 = []; %pk5 = []; pk6 = [];
lc = []; vr = [];
lc1 = []; lc3 = []; lc4 = []; lc2 = []; %lc5 = []; lc6 = [];
progressbar(0.3)
for i=1:k-1
    pr = detrend(pres(brkpoint(i):brkpoint(i+1)));
    pos=brkpoint(i)-1;
    pres_sf = sgolayfilt(pr,11,1+1000); % 1+1000 au lieu de 1+2*Fs ca marche mieux
    vol = detrend(cumsum(pres_sf))/1000;
    MPH = mean(abs(pres_sf))/2;
    MPH3 = mean(abs(vol))/2;
    
    if strcmp(mark,'peakmax')
        [pks1,locs1] = findpeaks(pres_sf,'MINPEAKHEIGHT',MPH,'MINPEAKDISTANCE',MINPEAKDISTANCE);
        % Variability signals
        vr1 = diff(locs1);
        ind1=[];
        for i1=6:length(vr1)
            if vr1(i1) < mean(vr1(i1-5:i1-1))*0.6 ||  0.2*vr1(i1) >mean(vr1(i1-5:i1-1))
                ind1=[ind1 i1];
            end;
        end;
        locs1(ind1)=[]; pks1(ind1)=[];
        pk1 = [pk1 pks1];
        lc1 = [lc1 locs1+pos];
        lc = lc1;
        vr = vr1;
        
        % Plot
%         for i=1:length(lc1)
%             vline(lc1(i),'r');
%         end
        
    elseif strcmp(mark,'peakmin')
        [pks2,locs2] = findpeaks(-pres_sf,'MINPEAKHEIGHT',MPH,'MINPEAKDISTANCE',MINPEAKDISTANCE);
        % Variability signals
        vr2 = diff(locs2);
        ind2=[];
        for i2=6:length(vr2)
            if vr2(i2) < mean(vr2(i2-5:i2-1))*0.6 ||  0.2*vr2(i2) >mean(vr2(i2-5:i2-1))
                ind2=[ind2 i2];
            end;
        end;
        locs2(ind2)=[]; pks2(ind2)=[];
        pk2 = [pk2 pks2];
        lc2 = [lc2 locs2+pos];
        lc = lc2;
        vr = vr2;
        
        % Plot
%         for i=1:length(lc2)
%             vline(lc2(i),'g');
%         end
        
        
    elseif strcmp(mark,'startexhal')
        [pks3,locs3] = findpeaks(vol,'MINPEAKDISTANCE',MINPEAKDISTANCE);
        % Variability signals
        vr3 = diff(locs3);
        ind3=[];
        for i3=6:length(vr3)
            if vr3(i3) < mean(vr3(i3-5:i3-1))*0.6 ||  0.2*vr3(i3) >mean(vr3(i3-5:i3-1))
                ind3=[ind3 i3];
            end
        end
        locs3(ind3)=[]; pks3(ind3)=[];
        pk3 = [pk3 pks3];
        lc3 = [lc3 locs3+pos];
        lc = lc3;
        vr = vr3;
        % Plot
        %figure;
        %plot(cumsum(detrend(pres))/100);
        
%         for i=1:length(lc3)
%             vline(lc3(i),'r');
%         end
        
        %         vol = detrend(cumsum(-pres_sf))/1000;
        %         MPH3 = mean(abs(vol))/2;
        %         [pks5,locs5] = findpeaks(vol,'MINPEAKDISTANCE',2*Fs);
        %         % Variability signals
        %         vr5 = diff(locs5);
        %         ind5=[];
        %         for i5=6:length(vr5)
        %             if vr5(i5) < mean(vr5(i5-5:i5-1))*0.6 ||  0.2*vr5(i5) >mean(vr5(i5-5:i5-1))
        %                 ind5=[ind5 i5];
        %             end
        %         end
        %         locs5(ind5)=[]; pks5(ind5)=[];
        %         pk5 = [pk5 pks5];
        %         lc5 = [lc5 locs5+pos];
        %         lc = lc5;
        %         % Plot
        %         %figure;
        %         %plot(cumsum(detrend(pres))/100);
        %         for i=1:length(lc5)
        %             vline(lc5(i),'g');
        %         end
        
    elseif strcmp(mark, 'startinspi')
        [pks4,locs4] = findpeaks(-vol,'MINPEAKDISTANCE',MINPEAKDISTANCE);
        % Variability signals
        vr4 = diff(locs4);
        ind4=[];
        for i4=6:length(vr4)
            if vr4(i4) < mean(vr4(i4-5:i4-1))*0.6 ||  0.2*vr4(i4) >mean(vr4(i4-5:i4-1))
                ind4=[ind4 i4];
            end
        end
        locs4(ind4)=[]; pks4(ind4)=[];
        pk4 = [pk4 pks4];
        lc4 = [lc4 locs4+pos];
        lc = lc4;
        vr = vr4;
        % Plot
        %         figure;
        %         plot(cumsum(detrend(pres))/100);
        
%         for i=1:length(lc4)
%             vline(lc4(i),'r');
%         end;
        
        %         vol = detrend(cumsum(-pres_sf))/1000;
        %         MPH3 = mean(abs(vol))/2;
        %         [pks6,locs6] = findpeaks(-vol,'MINPEAKDISTANCE',2*Fs);
        %         % Variability signals
        %         vr6 = diff(locs6);
        %         ind6=[];
        %         for i6=6:length(vr6)
        %             if vr6(i6) < mean(vr6(i6-5:i6-1))*0.6 ||  0.2*vr6(i6) >mean(vr6(i6-5:i6-1))
        %                 ind6=[ind6 i6];
        %             end
        %         end
        %         locs6(ind6)=[]; pks6(ind6)=[];
        %         pk6 = [pk6 pks6];
        %         lc6 = [lc6 locs6+pos];
        %         lc = lc6;
        %         % Plot
        %         %figure;
        %         %plot(cumsum(detrend(pres))/100);
        %         for i=1:length(lc6)
        %             vline(lc6(i),'g');
        %         end
    end;
	if k-1==1
		progressbar(0.45)
		progressbar(0.6)
		progressbar(0.75)
		progressbar(0.9)
	 else
		if i==floor((k-1)/6)
			progressbar(0.4)
		elseif i==floor(2*(k-1)/6)
			progressbar(0.5)
		elseif i==floor(3*(k-1)/6)
			progressbar(0.6)
		elseif i==floor(4*(k-1)/6)
			progressbar(0.7)
		elseif i==floor(5*(k-1)/6)
			progressbar(0.8)
		elseif i==k-1
			progressbar(0.9)
		end;
	 end;
    %disp(['j = ' num2str(i) ])
end
progressbar(1)

disp('Done.');
return
