function distance = featureBasedDTW(ref, test, slope)

% Implements Feature Based Dynamic Time Warping
% ref: reference MFCC vector with dimensions numCoefficients by numFrames
% test: test MFCC vector with dimensions numCoefficients by numFrames
% slope: slope for the parallelogram 

if nargin<3; slope=2; end

% Instead of taking the raw MFCC coefficients, we estimate the derivative
% to obtain more refined features

[ref_numCoef ref_numFrames] = size(ref);
[test_numCoef test_numFrames] = size(test);

% Local and Global features for ref vector
ref_localX = zeros(ref_numCoef, ref_numFrames-2);
ref_localY = zeros(ref_numCoef, ref_numFrames-2);
ref_globalX = zeros(ref_numCoef, ref_numFrames-2);
ref_globalY = zeros(ref_numCoef, ref_numFrames-2);
for i = 1:ref_numCoef
    for j = 2:ref_numFrames-1
        ref_localX(i, j-1) = ref(i, j) - ref(i, j-1);
        ref_localY(i, j-1) = ref(i, j) - ref(i, j+1);
        sumX = 0;
        for k = 1:j-1
            sumX = sumX + ref(i, k);
        end
        sumX = sumX/(j-1);
        ref_globalX(i, j-1) = ref(i, j) - sumX;
        sumY = 0;
        for k = j+1:ref_numFrames
            sumY = sumY + ref(i, k);
        end
        sumY = sumY/(ref_numFrames-j);
        ref_globalY(i, j-1) = ref(i, j) - sumY;
    end
end

% Local and Global features for test vector
test_localX = zeros(test_numCoef, test_numFrames-2);
test_localY = zeros(test_numCoef, test_numFrames-2);
test_globalX = zeros(test_numCoef, test_numFrames-2);
test_globalY = zeros(test_numCoef, test_numFrames-2);
for i = 1:test_numCoef
    for j = 2:test_numFrames-1
        test_localX(i, j-1) = test(i, j) - test(i, j-1);
        test_localY(i, j-1) = test(i, j) - test(i, j+1);
        sumX = 0;
        for k = 1:j-1
            sumX = sumX + test(i, k);
        end
        sumX = sumX/(j-1);
        test_globalX(i, j-1) = test(i, j) - sumX;
        sumY = 0;
        for k = j+1:test_numFrames
            sumY = sumY + test(i, k);
        end
        sumY = sumY/(test_numFrames-j);
        test_globalY(i, j-1) = test(i, j) - sumY;
    end
end

ref_numFrames = ref_numFrames - 2;
test_numFrames = test_numFrames - 2;

% Sakoe - Chiba band
R = max(ref_numFrames, test_numFrames)*.25; 
R = ceil(R);

% Identify the region for the optimum warping path
% Slop decides the parallelogram dimensions; for Itakura use 2 & 1/2; one
% paper suggests using 3 & 1/3
I = ref_numFrames; J = test_numFrames;
flag = zeros(I, J);
for i = 1:ref_numFrames
    for j = 1:test_numFrames
        % Sakoe-Chiba band
        if (abs(j - (J/I)*i) <= R) 
            flag(i, j) = 1;
        end
        % Slope band
        if ((abs(j - slope*i) <= abs(J - slope*I)) && (abs(j - i/slope) <= abs(J - I/slope)))
            flag(i, j) = 1;
        end
    end
end

% Compute the distance between the MFCC vectors in this region
local_dist = zeros(ref_numFrames, test_numFrames);

% Initialize to Infinity (and not to zero, since below we will be taking the minimum
% distance)
for i=1:ref_numFrames
    for j=1:test_numFrames     
        local_dist(i, j) = inf;
    end
end
sum = 0;

for i=1:ref_numFrames
    for j=1:test_numFrames    
        if (flag(i, j) == 1) % Compute only for the valid region
            for k=1:ref_numCoef
                distlocal = abs(ref_localX(k,i) - test_localX(k,j)) + abs(ref_localY(k,i) - test_localY(k,j));
                distglobal = abs(ref_globalX(k,i) - test_globalX(k,j)) + abs(ref_globalY(k,i) - test_globalY(k,j));
                sum = sum + distlocal + distglobal;
            end
            local_dist(i,j) = sum;
        end
        sum = 0;
    end
end

accum_dist(1,1) = local_dist(1,1);

    for x = 2:test_numFrames
    
        accum_dist(1,x) = accum_dist(1,x-1)+local_dist(1,x);
    
    end
    
    for y = 2:ref_numFrames
    
        accum_dist(y,1) = accum_dist(y-1,1)+local_dist(y,1);
    
    end
    
for a = 2:ref_numFrames
    for b = 2:test_numFrames
        
        accum_dist(a,b) = local_dist(a,b) + min([accum_dist(a-1,b) accum_dist(a,b-1) accum_dist(a-1,b-1)]);
        
    end
end

%distance = accum_dist(ref_numFrames, test_numFrames);
distance = accum_dist(ref_numFrames, test_numFrames)/(ref_numFrames + test_numFrames);

