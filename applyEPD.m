% Format of training filename: numSpeaker_textDigit_numEmotion.wav
% eg. 7_text2_6.wav => represents speaker number 7, speaking text 2, in
% Neutral emotion

% Codes for Emotion
% 1: Anger; 2: Boredom; 3: Fear/Anxiety; 4: Happiness; 5: Sadness; 6: Neutral

% Number of Speakers from 7 to 16 (10 speakers)

clc;clear;close all;

for numSpeaker = 7:16
    for textCode = 1:10
        for numEmotion = 1:6
            reffile = sprintf('Database/%d_text%d_%d.wav', numSpeaker, textCode, numEmotion);
%             if (reffile=='Database/8_text6_1.wav')
%                 a=1;
%             end
            [refinput fs] = wavread(reffile);
            [ITL,ITU,IZCT,energy,zc,N1,N2]=epdStats(refinput, fs);
            if (N2 > length(refinput))
                N2=length(refinput);
            end
            refinput = refinput(N1:N2);
            outfile = sprintf('Database/%d_text%d_%d_epd.wav', numSpeaker, textCode, numEmotion);
            wavwrite(refinput, fs, 16, outfile);
        end
    end
end
