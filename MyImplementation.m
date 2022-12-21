% Format of training filename: numSpeaker_digit_numUtterance
% eg. 1_1_1.wav, 1_1_2.wav, 1_1_3.wav, 1_2_1.wav, 1_4_3.wav, 3_1_1.wav etc.
% Format of test filename: test_numSpeaker_digit_numUtterance

clc;clear;close all;

counter = 0;
accuracyCount = 0;
totalNumUtterances = 12;
spokenDigit = 1; % 0 to 9

for spokenDigit = 0:9
    t0 = clock;
    
    for testDigit = 1:8
        distance = zeros(totalNumUtterances*10, 1);
        
        testfile = sprintf('TIDataEPD/Test/test_%d_%d_epd.wav', spokenDigit, testDigit);   
        [testinput fs1] = wavread(testfile);
        test_mfcc = wave2mfcc(testinput);
        
        for digit = 0:9
            for numUtterance = 1:totalNumUtterances
                counter = counter + 1;
                reffile = sprintf('TIDataEPD/Train/%d_%d_epd.wav', digit, numUtterance);
                [refinput fs2] = wavread(reffile);
                ref_mfcc = wave2mfcc(refinput);
                distance(counter, 1) = myDTW(test_mfcc, ref_mfcc); 
            end
        end
        
        [minimum index] = min(distance);
        match  = floor((index-1)/double(totalNumUtterances));
        if (match == spokenDigit)
            accuracyCount = accuracyCount + 1;
        end
        sprintf('For test digit %d_%d, the detected number is %d    ', spokenDigit, testDigit, match)
        
        
        counter = 0;
    end
    
    %sprintf('Time elapsed for digit %d is %f', spokenDigit, etime(clock, t0))
    
end

sprintf('Overall Accuracy = %f',(accuracyCount/80)*100)


