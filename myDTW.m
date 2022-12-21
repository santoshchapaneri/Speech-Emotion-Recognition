function [distance, dtwPath] = myDTW(ref, test, typeDTW, slope)

% typeDTW = Type of DTW 0: DTW, 1: DDTW, 2: FBDTW, 3: IFDTW

% Parameter Set for DTW
DP = dtwPrmSet;

if nargin<3
    typeDTW = DP.typeDTW;
end

if nargin<4 
    slope = DP.slope; 
end

[ref_numCoef ref_numFrames] = size(ref);
[test_numCoef test_numFrames] = size(test);

if (typeDTW == 1) % DDTW
    % Instead of taking the raw MFCC coefficients, we estimate the derivative
    % to obtain more refined features
    [ref, test] = computeDerivativeEstimates(ref, test);
end

if (typeDTW == 2) % FBDTW
    [ref_localX, ref_localY, ref_globalX, ref_globalY, test_localX, test_localY, test_globalX, test_globalY] = computeFBDTW(ref, test);
    ref_numFrames = ref_numFrames - 2;
    test_numFrames = test_numFrames - 2;
end

if (typeDTW == 3) % IFDTW
    [ref_local, ref_global, test_local, test_global] = computeIFDTW(ref, test);
    ref_numFrames = ref_numFrames - 2;
    test_numFrames = test_numFrames - 2;
end

I = ref_numFrames; J = test_numFrames;
flag = ones(I, J); % region for parallelogram - default is full region of I x J

% Find the region for warping path - global constraint
if (DP.SCBand | DP.slopeBand)
    flag = zeros(I, J); % region for parallelogram - init to 0
    
    % Sakoe - Chiba band
    if (DP.SCBand)
        width = DP.SCWidth;
        R = max(ref_numFrames, test_numFrames)*width; 
        R = ceil(R);
    end

    % Identify the region for the optimum warping path
    % Slope decides the parallelogram dimensions; for Itakura use 2 & 1/2; one
    % paper suggests using 3 & 1/3
    for i = 1:I
        for j = 1:J
            % Sakoe-Chiba band
            if (DP.SCBand)
                if (abs(j - (J/I)*i) <= R) 
                    flag(i, j) = 1;
                end
            end
            % Slope band
            if(DP.slopeBand)
                if ((i - j/slope >= 0) && (i - slope*j <= 0) && (i - slope*j >= I - slope*J) && (i - j/slope <= I - J/slope))
                    flag(i, j) = 1;
                end
            end
        end
    end
end


% Compute the distance between the input MFCC vectors
local_dist = zeros(ref_numFrames, test_numFrames);

% Initialize to Infinity (and not to zero, as we will be taking the minimum distance)
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
                
                if (typeDTW == 0 | typeDTW == 1) % Regular DTW or Derivative DTW
                    if (DP.metric == 1)
                        sum = sum + (abs(ref(k,i)-test(k,j))); % City-block distance
                    else 
                        sum = sum + (ref(k,i)-test(k,j))^2; % Euclidean distance
                    end
                end
                
                if (typeDTW == 2) % Feature Based DTW
                    distlocal = abs(ref_localX(k,i) - test_localX(k,j)) + abs(ref_localY(k,i) - test_localY(k,j));
                    distglobal = abs(ref_globalX(k,i) - test_globalX(k,j)) + abs(ref_globalY(k,i) - test_globalY(k,j));
                    sum = sum + distlocal + distglobal;
                end
                
                if (typeDTW == 3) % Improved Features DTW
                    distlocal = abs(ref_local(k,i) - test_local(k,j));
                    distglobal = abs(ref_global(k,i) - test_global(k,j));
                    sum = sum + distlocal + distglobal;
                end
                
            end
            if (DP.metric == 1 | typeDTW == 2 | typeDTW == 3)
                local_dist(i,j) = sum;
            else
                local_dist(i,j) = sqrt(sum); % Euclidean Distance
            end
        end
        sum = 0;
    end
end

% ====== Construct prevPos table for back tracking the optimum path
for i=1:ref_numFrames
    for j=1:test_numFrames
		prevPos(i,j).i=-1;
		prevPos(i,j).j=-1;
	end
end

% Accumulated distance over the optimum warping path
accum_dist(1,1) = local_dist(1,1);

for x = 2:test_numFrames
    accum_dist(1, x) = accum_dist(1, x-1)+local_dist(1, x);
end
    
for y = 2:ref_numFrames
    accum_dist(y, 1) = accum_dist(y-1, 1)+local_dist(y, 1);
end
    
for i = 2:ref_numFrames
    for j = 2:test_numFrames
        minValue = min([accum_dist(i-1, j) accum_dist(i, j-1) accum_dist(i-1, j-1)]);
        accum_dist(i, j) = local_dist(i, j) + minValue;

        if (minValue == accum_dist(i-1, j-1))
            prevPos(i,j).i=i-1; prevPos(i,j).j=j-1; 
        else if (minValue == accum_dist(i-1, j))
            prevPos(i,j).i=i-1; prevPos(i,j).j=j; 
        else
            prevPos(i,j).i=i; prevPos(i,j).j=j-1; 
        end
        
    end
end

% Final distance (normalized)
distance = accum_dist(ref_numFrames, test_numFrames)/(ref_numFrames + test_numFrames);

% ====== Back track to find the warping path
dtwPath=[ref_numFrames; test_numFrames];		% The last point in the optimum path
nextPoint=[prevPos(dtwPath(1,1), dtwPath(2,1)).i; prevPos(dtwPath(1,1), dtwPath(2,1)).j];
while nextPoint(1)>0 & nextPoint(2)>0
	dtwPath=[nextPoint, dtwPath];
	nextPoint=[prevPos(dtwPath(1,1), dtwPath(2,1)).i; prevPos(dtwPath(1,1), dtwPath(2,1)).j];
end

end % of main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes derivative estimates for Derivative Dynamic Time Warping
function [ref_derivative, test_derivative] = computeDerivativeEstimates(ref, test)
[ref_numCoef ref_numFrames] = size(ref);
[test_numCoef test_numFrames] = size(test);

% Derivative estimates for ref vector
ref_derivative = zeros(ref_numCoef, ref_numFrames);
for i = 1:ref_numCoef
    for j = 2:ref_numFrames-1
        ref_derivative(i, j) = ( (ref(i, j) - ref(i, j-1)) + ((ref(i, j+1) - ref(i, j-1))/2) ) / 2;
    end
end
% Copy over the second and penultimate elements to first and last respectively
ref_derivative(:,1) = ref_derivative(:,2);
ref_derivative(:,ref_numFrames) = ref_derivative(:,ref_numFrames-1);

% Derivative estimates for test vector
test_derivative = zeros(test_numCoef, test_numFrames);
for i = 1:test_numCoef
    for j = 2:test_numFrames-1
        test_derivative(i, j) = ( (test(i, j) - test(i, j-1)) + ((test(i, j+1) - test(i, j-1))/2) ) / 2;
    end
end
% Copy over the second and penultimate elements to first and last respectively
test_derivative(:,1) = test_derivative(:,2);
test_derivative(:,test_numFrames) = test_derivative(:,test_numFrames-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements Feature Based Dynamic Time Warping
function [ref_localX, ref_localY, ref_globalX, ref_globalY, test_localX, test_localY, test_globalX, test_globalY] = computeFBDTW(ref, test)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements Improved Features Dynamic Time Warping (IFDTW)
function [ref_local, ref_global, test_local, test_global] = computeIFDTW(ref, test)
[ref_numCoef ref_numFrames] = size(ref);
[test_numCoef test_numFrames] = size(test);

% Local and Global features for ref vector
ref_local = zeros(ref_numCoef, ref_numFrames-2);
ref_global = zeros(ref_numCoef, ref_numFrames-2);
for i = 1:ref_numCoef
    for j = 2:ref_numFrames-1
        ref_local(i, j-1) = ref(i, j) - (ref(i, j+1) - ref(i, j-1))/2;
        sum1 = 0;
        for k = j+1:ref_numFrames
            sum1 = sum1 + ref(i, k);
        end
        sum1 = sum1/(ref_numFrames-j);

        sum2 = 0;
        for k = 1:j-1
            sum2 = sum2 + ref(i, k);
        end
        sum2 = sum2/(j-1);
        ref_global(i, j-1) = ref(i, j) - (sum1 - sum2)/2;
    end
end

% Local and Global features for test vector
test_local = zeros(test_numCoef, test_numFrames-2);
test_global = zeros(test_numCoef, test_numFrames-2);
for i = 1:test_numCoef
    for j = 2:test_numFrames-1
        test_local(i, j-1) = test(i, j) - (test(i, j+1) - test(i, j-1))/2;
        sum1 = 0;
        for k = j+1:test_numFrames
            sum1 = sum1 + test(i, k);
        end
        sum1 = sum1/(test_numFrames-j);

        sum2 = 0;
        for k = 1:j-1
            sum2 = sum2 + test(i, k);
        end
        sum2 = sum2/(j-1);
        test_global(i, j-1) = test(i, j) - (sum1 - sum2)/2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
