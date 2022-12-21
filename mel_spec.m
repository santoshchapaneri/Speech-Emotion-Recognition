function [melspec] = mel_spec(x, Nfft, winl, hp, gH)

% x:  input speech data
% Nfft: No. of FFT points
% winl:  analysis window length
% hp:  frame shift
% gH : Matrix of Mel-Filter weights

% hamming window based analysis
hwin = hamming(winl);

x = filter([1 -0.97], 1, x);

f = aSTFT(x, Nfft, winl, hp, hwin);
spec = abs(f).^2; 
melspec = gH * spec;


function [P] = aSTFT(x, Nfft, winl, hp, hwin)
x = x(:);
% L = winl + (N-1)*hp
N = ceil((length(x) - winl)/hp + 1);
Nzeros = winl + (N-1)*hp - length(x);

x = [x; zeros(Nzeros, 1)];

s = zeros(winl, N);
i = winl;
for j = 1:N
    s(:,j) = hwin .* x((i-winl+1):i);
    i = i + hp;
end

% FFT it
f = fft(s, Nfft);

% Chop redundant part
P = f(1:end/2+1,:);
P = abs(P);


function [P] = sSTFT(P)
P = P;