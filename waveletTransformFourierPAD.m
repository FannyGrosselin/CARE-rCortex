function [Z,freqVector] = waveletTransformFourierPAD(f1,f2,pas,nfft,X)
%
% wavelet transform (Morlet complex wavelet)  of signal X
%
% Z = waveletTransformFourier(f1,f2,pas,k,x)
%
% Author: Mario Chavez
%
%  INPUTS:
%      X         is a column  vector
%      f1        is the lower frequency in the decomposition (in normalized units, ie the highest frequency is 0.5)
%      f2        is the higher frequency in the decomposition (in normalized units, ie the highest frequency is 0.5)
%      pas      frequency resolution (in normalized units)
%
%  OUTPUT
%      Z       wavelet transform of X (can be complex valued)

n = length(X);
X = X - mean(X);
if n == size(X,2)
    X = X';
end;
if nfft>n
    X = [X; zeros(nfft-n,1)];
end
freqVector = [f1:pas:f2]';

widthMorlet = 5;   %--- width of the mother wavelet
m = length(f1:pas:f2);
l=1;

taille = max(size(X));
if rem(taille,2)
    Nfref = (taille+1)/2;
else
    Nfref = taille/2+1;
end

freq = ((0:(Nfref-1))/(Nfref-1))/2;
fftx = fft(X);
Psi = zeros(1,taille);
Z = zeros(m,n);
for f = f1:pas:f2
    sigmaF = f/widthMorlet;
	if sigmaF ~= 0
		sigmaT = 1/(2*pi*sigmaF);
	else
		sigmaT = 0;
	end;
    w = 2*pi*(freq - f)*sigmaT;
    Psi(1:Nfref) = realpow(4*pi*sigmaT*sigmaT, 1/4)*exp(-(w .* w)/2);
    %figure, plot(Psi), pause;
    dummy = fliplr(ifft(fftx'.* Psi));
    Z(l,:)  = dummy(1:n);
    l = l+1;
end
return