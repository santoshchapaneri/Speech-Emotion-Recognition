function distance = derivativeDTW(ref, test, slope)

% Implements Derivative Dynamic Time Warping
% ref: reference MFCC vector with dimensions numCoefficients by numFrames
% test: test MFCC vector with dimensions numCoefficients by numFrames
% slope: slope for the parallelogram 

if nargin<3; slope=2; end

% Instead of taking the raw MFCC coefficients, we estimate the derivative
% to obtain more refined features

[ref_numCoef ref_numFrames] = size(ref);
[test_numCoef test_numFrames] = size(test);

% Derivative estimates for ref vector
ref_derivative = zeros(ref_numCoef, ref_numFrames);
for i = 1:ref_numCoef
    for j = 2:ref_numFrames-1
        ref_derivative(i, j) = ( (ref(i, j) - ref(i, j-1)) + ((ref(i, j+1) - ref(i, j-1))/2) ) / 2;
    end
end
% Copy over the second and penultimate elements to first and last
% respectively
ref_derivative(:,1) = ref_derivative(:,2);
ref_derivative(:,ref_numFrames) = ref_derivative(:,ref_numFrames-1);

% Derivative estimates for test vector
test_derivative = zeros(test_numCoef, test_numFrames);
for i = 1:test_numCoef
    for j = 2:test_numFrames-1
        test_derivative(i, j) = ( (test(i, j) - test(i, j-1)) + ((test(i, j+1) - test(i, j-1))/2) ) / 2;
    end
end
% Copy over the second and penultimate elements to first and last
% respectively
test_derivative(:,1) = test_derivative(:,2);
test_derivative(:,test_numFrames) = test_derivative(:,test_numFrames-1);

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
                sum = sum + (abs(ref_derivative(k,i)-test_derivative(k,j))); % City-block distance
                %sum = sum + (t(i,k)-r(j,k))^2; % Euclidean distance
            end
        %local_dist(i,j) = sqrt(sum); %For Euclidean Distance
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

