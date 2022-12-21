vec1=[71 73 75 80 80 80 78 76 75 73 71 71 71 73 75 76 76 68 76 76 75 73 71 70 70 69 68 68 72 74 78 79 80 80 78];
vec2=[69 69 73 75 79 80 79 78 76 73 72 71 70 70 69 69 69 71 73 75 76 76 76 76 76 75 73 71 70 70 71 73 75 80 80 80 78];

train_4_1 = wavread('TIDataEPD\Train\4_1_epd.wav');
train_mfcc_4_1 = wave2mfcc(train_4_1);
train_5_1 = wavread('TIDataEPD\Train\5_1_epd.wav');
train_mfcc_5_1 = wave2mfcc(train_5_1);
train_5_3 = wavread('TIDataEPD\Train\5_3_epd.wav');
train_mfcc_5_3 = wave2mfcc(train_5_3);
train_7_8 = wavread('TIDataEPD\Train\7_8_epd.wav');
train_mfcc_7_8 = wave2mfcc(train_7_8);
test_4_2 = wavread('TIDataEPD\Test\test_4_2_epd.wav');
test_mfcc_4_2 = wave2mfcc(test_4_2);
test_4_8 = wavread('TIDataEPD\Test\test_4_8_epd.wav');
test_mfcc_4_8 = wave2mfcc(test_4_8);

vec1 = train_mfcc_4_1(:,6)';
vec2 = test_mfcc_4_8(:,6)';

slope = 2;
[minDist0, dtwPath0] = myDTW(vec1, vec2, 0, slope);
[minDist1, dtwPath1] = myDTW(vec1, vec2, 1, slope);
[minDist2, dtwPath2] = myDTW(vec1, vec2, 2, slope);
[minDist3, dtwPath3] = myDTW(vec1, vec2, 3, slope);
% subplot(4,1,1); dtwBridgePlot(vec1, vec2, dtwPath0, '2d'); title('Alignment by regular DTW');
% subplot(4,1,2); dtwBridgePlot(vec1, vec2, dtwPath1, '2d'); title('Alignment by derivative DTW');
% subplot(4,1,3); dtwBridgePlot(vec1, vec2, dtwPath2, '2d'); title('Alignment by FBDTW');
% subplot(4,1,4); dtwBridgePlot(vec1, vec2, dtwPath3, '2d'); title('Alignment by IFDTW');
subplot(3,1,1); dtwBridgePlot(vec1, vec2, dtwPath0, '2d'); title('Alignment by DTW');
subplot(3,1,2); dtwBridgePlot(vec1, vec2, dtwPath1, '2d'); title('Alignment by Derivative DTW');
subplot(3,1,3); dtwBridgePlot(vec1, vec2, dtwPath3, '2d'); title('Alignment by Improved Features DTW');

% figure; dtwPathPlot(vec1, vec2, dtwPath0); title('Alignment by regular DTW');
% figure; dtwPathPlot(vec1, vec2, dtwPath1); title('Alignment by DDTW');
% figure; dtwPathPlot(vec1, vec2, dtwPath2); title('Alignment by FBDTW');
% figure; dtwPathPlot(vec1, vec2, dtwPath3); title('Alignment by IFDTW');



