function [] = features_MFCC_DSCC(wavfile, featfile)

% Code for speech feature (MFCC and DSCC) extraction

% wavfile: speech wav file in MS-WAV format
% featfile: file to write obtained features in Sphinx feature format

% feature parameters
fs = 16000;        % Sampling frequency
winl = 0.025 * fs; % Window length (25 ms)
hp   = 0.010 * fs; % Frame shift (10 ms)
Nfilt = 40;        % No. of filters in Mel-Filter bank
Nfft = 512;        % No. of FFT points
Ncep = 13;         % No. of Cepstral Coeffs.
T = 3;             % parameter for DSCC

% Assign the temppath variable, will be used for creating intermediate
% files
temppath = './';

melfile = sprintf('%s/melfilt_Nfilt%dNfft%d.mat', temppath, Nfilt, Nfft);
if(exist(melfile,'file') & 0);
    load(melfile,'gH');
else
    gH = get_melfilt(Nfilt, Nfft+2);
    save(melfile,'gH');
end

% read wav file
input = wavread(wavfile);

% data normalization
input = input - mean(input);     input = input/std(input);

% Mel spectral features
P = mel_spec(input, Nfft, winl, hp, gH);

nonlP = log(P);
% MFCC - Keep Ncep no. of DCT features
MFCC = dct(nonlP);
MFCC = MFCC(1:Ncep,:);
MFCC = mean_var_norm(MFCC);

% DSCC - included delta-spec and double-delta-spec features
DSCC = specD(P, T, Ncep);

feats = [MFCC; DSCC];
feats = feats(:,2:end-1); % Ignoring the first and last frame [due to DSCC processing]
writemfcc(featfile, feats);