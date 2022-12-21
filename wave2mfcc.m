function [parameter]= wave2mfcc(input, fs, FP)
% wave2mfcc: Wave to MFCC (Mel-Frequency Cepstral Cofficient) conversion
%	Usage:
%	parameter = wave2mfcc(input, fs, FP)
%	parameter: MFCC and log energy, plus their delta value if necessary.
%   input: speech signal (user must prior call wavread on the wav file)
%	fs: sampling rate
%	FP: Parameters for deriving the speech features. You can use mfccPrmSet.m to obtain the parameters.
%
%	For example:
%		waveFile='speechfilename.wav';
%		[input, fs]=wavread(waveFile);
%		FP=mfccPrmSet(fs);
%		FP.useDelta=0;
%		mfcc0= wave2mfcc(input, fs, FP);
%		fprintf('No. of extracted frames = %d\n', size(mfcc0, 2));
%		subplot(3,1,1); surf(mfcc0); box on; axis tight; title(strPurify(sprintf('13D MFCC of "%s"', waveFile)));
%		FP.useDelta=1;
%		mfcc1=wave2mfcc(input, fs, FP);
%		subplot(3,1,2); surf(mfcc1); box on; axis tight; title(strPurify(sprintf('26D MFCC of "%s"', waveFile)));
%		FP.useDelta=2;
%		mfcc2=wave2mfcc(input, fs, FP);
%		subplot(3,1,3); surf(mfcc2); box on; axis tight; title(strPurify(sprintf('39D MFCC of "%s"', waveFile)));

if nargin<1, selfdemo; return; end
if nargin<2; fs=16000; end
%fs = 16000;
if nargin<3, FP=mfccPrmSet(fs); end

input=double(input);	% Convert to double
% input=input-mean(input);	% Shift to zero mean
% input = input./std(input);  % Make Std. deviation = 1

% ====== Step 1: Pre-Emphasis
inputPreEmp = filter([1, -FP.preEmCoef], 1, input);
%inputPreEmp = input;

% ====== Step 2: Frame Blocking
framedInput = buffer2(inputPreEmp, FP.frameSize, FP.overlap);

filterBankParam = getTriFilterPrm(FP.frameSize, fs, FP.tbfNum, 0);	% Parameters for triangular filter bank

parameter = [];
for i = 1:size(framedInput, 2),
	
    % ====== Step 3: Hamming Window
	Wframe  = hamming(FP.frameSize).*framedInput(:,i);
	
    % ====== Step 4: FFT
	fftMag = abs(fft(Wframe, FP.NumFft));
	halfIndex = floor((FP.NumFft+1)/2); % Spectrum is symmetric
	fftMag = fftMag(1:halfIndex);
	
    % ====== Step 5: Mel Triangular Bandpass Filter
	tbfCoef = triBandFilter(fftMag, FP.tbfNum, filterBankParam); % See below for function

    % ====== Step 6: DCT to get L order mel frequency cepstrum coefficients
	mfcc = melCepstrum(FP.cepsNum, FP.tbfNum, tbfCoef);
	
    % ====== Step 7: Liftering
    if (FP.UseLiftering)
        lifter=1:FP.cepsNum;    %Lifter vector index
        lifter=1+floor((FP.cepsNum)/2)*(sin(lifter*pi/FP.cepsNum));    %raised sine lifter version
        mfcc = mfcc.*lifter;    %Multiplies mfcc by lifter values
    end
	
    parameter = [parameter mfcc'];
    
end

% ====== Add energy coefficient
if (FP.useEnergy==1)
	energy = sum(framedInput.^2)/FP.frameSize;
	logEnergy = 10*log10(eps+energy);
	parameter = [parameter; logEnergy];
end

% Perform Cepstral Mean Substraction & Variance Normalization
if (FP.useCMS)
    % Mean normalization
    meanparam = mean(parameter, 2);
    parameter = parameter - repmat(meanparam, 1, size(parameter, 2));
    
    % Variance normalization
    M = size(parameter, 1);
    std_parm = ones(M, 1);
    vp = sqrt(mean(parameter .* parameter, 2));
    vp = vp./std_parm;
    parameter = parameter./repmat(vp, 1, size(parameter, 2));
end

% ====== Compute Delta Energy and Delta Cepstrum Coefficients
if (FP.useDelta>=1)
	deltaWindow = 2;
	paraDelta = deltaFunction(deltaWindow, parameter);
	parameter = [parameter; paraDelta];
end
if (FP.useDelta>=2)
	paraDeltaDelta = deltaFunction(deltaWindow, paraDelta);
	parameter = [parameter; paraDeltaDelta];
end
% Compute Triple Delta Features
if (FP.useDelta==3)
	paraDeltaDeltaDelta = deltaFunction(deltaWindow, paraDeltaDelta);
	parameter = [parameter; paraDeltaDeltaDelta];
end

% parameter has dimensions of NumCoeffs by numFrames

% Modified MFCC coefficients
if (FP.newMFCC && (FP.useDelta==1 | FP.useDelta==2 | FP.useDelta==3))
    a = 1/2;
    b = 1/4;
    c = 1/8;
    numMFCC = FP.cepsNum + FP.useEnergy;
    newMFCC = parameter(1:numMFCC,:) + a*parameter(numMFCC+1:2*numMFCC,:) + b*parameter(2*numMFCC+1:3*numMFCC,:) + c*parameter(3*numMFCC+1:4*numMFCC,:);
    parameter = newMFCC;
end

% ====== Subfunction ======
% === Self demo
function selfdemo
waveFile='Data\0_1.wav';
[input, fs]=wavread(waveFile);
fs=16000;
FP=mfccPrmSet(fs);
FP.useDelta=0;
mfcc0= feval(mfilename, input, fs, FP);
fprintf('No. of extracted frames = %d\n', size(mfcc0, 2));
subplot(3,1,1); surf(mfcc0); box on; axis tight; title(strPurify(sprintf('13D MFCC of "%s"', waveFile)));
FP.useDelta=1;
mfcc1=feval(mfilename, input, fs, FP);
subplot(3,1,2); surf(mfcc1); box on; axis tight; title(strPurify(sprintf('26D MFCC of "%s"', waveFile)));
FP.useDelta=2;
mfcc2=feval(mfilename, input, fs, FP);
subplot(3,1,3); surf(mfcc2); box on; axis tight; title(strPurify(sprintf('39D MFCC of "%s"', waveFile)));

% === Triangular band-pass filters
function tbfCoef = triBandFilter(fftMag, P, filterBankParam)
fstart=filterBankParam(1,:);
fcenter=filterBankParam(2,:);
fstop=filterBankParam(3,:);
% Triangular bandpass filter.
for i=1:P
   for j = fstart(i):fcenter(i),
      filtmag(j) = (j-fstart(i))/(fcenter(i)-fstart(i));
   end
   for j = fcenter(i)+1:fstop(i),
      filtmag(j) = 1-(j-fcenter(i))/(fstop(i)-fcenter(i));
   end
   tbfCoef(i) = sum(fftMag(fstart(i):fstop(i)).*filtmag(fstart(i):fstop(i))');
end
tbfCoef=log(eps+tbfCoef.^2); % eps to avoid log(0)

% === TBF coefficients to MFCC
function mfcc = melCepstrum(L, P, tbfCoef)
% DCT to find MFCC
for i = 1:L
   coef = cos((pi/P)*i*(linspace(1,P,P)-0.5))';
   mfcc(i) = sum(coef.*tbfCoef');
end

% === Delta function
function parameter = deltaFunction(deltaWindow, parameter)
% Compute Delta Cepstrum coefficients and Delta Log Energy
rows  = size(parameter,1);
cols  = size(parameter,2);
%temp  = [zeros(rows,deltaWindow) parameter zeros(rows,deltaWindow)];
temp  = [parameter(:,1)*ones(1,deltaWindow) parameter parameter(:,end)*ones(1,deltaWindow)];
temp2 = zeros(rows,cols);
denominator = sum([1:deltaWindow].^2)*2;
for i = 1+deltaWindow : cols+deltaWindow,
   subtrahend = 0;
   minuend    = 0;
   for j = 1 : deltaWindow,
      subtrahend = subtrahend + temp(:,i+j)*j;
      minuend = minuend + temp(:,i-j)*(-j);
   end;
   temp2(:,i-deltaWindow) = (subtrahend + minuend)/denominator;
end;
parameter = temp2;