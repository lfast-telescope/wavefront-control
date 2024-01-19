function [closestPoint, ptIndex] = FindClosestPoint(pointSet, point1, point2, point3)

% This function is used to find the point closest to the reference points.
% INPUTS:
%   pointSet: N x 2 point set in format [x y] used to find the point
% closest to reference points;
%   point1/point2/point3: reference points.
% OUTPUT: 
%   closestPoint: closest point;
%   ptIndex: closest point index in 'pointSet'.
% HISTORY:
%   2023-08-11 - Yiyang Huang - initial implementation

% Validate inputs.
if nargin < 2, error('Please give at least 1 point!'); end
if isempty(pointSet)
    closestPoint = nan; ptIndex = nan;
    return
end

% Find the closest point.
if nargin < 3 % 1 point
    ptIndex = dsearchn(pointSet,point1);
    closestPoint = pointSet(ptIndex,:);
elseif nargin < 4 % 2 points
    % Find the point which has the minimum sum of distances to point 1 & 2.
    ptsNum = size(pointSet,1); % unlabelled points number in the selected region
    centroidsTemp = pointSet; % temporary matrix used in the following processing
    k1 = []; k2 = [];
    for i = 1:6 % find 6 closest points to point 1 and 2 respectively
        if ptsNum == 0 % no point available, so stop the search
            break;
        end
        k1(i) = dsearchn(centroidsTemp,point1);
        k2(i) = dsearchn(centroidsTemp,point2);
        centroidsTemp(k1(i),:) = -100; centroidsTemp(k2(i),:) = -100; % set selected points to be invalid
        if k1(i) == k2(i)
            ptsNum = ptsNum - 1;
        else
            ptsNum = ptsNum - 2;
        end
    end
    k = union(k1,k2); % indices of points to be judged
    ptsNumTBD = length(k);
    for i = 1:ptsNumTBD % calculate sums of distances
        distance1(i) = norm(pointSet(k(i),:)-point1);
        distance2(i) = norm(pointSet(k(i),:)-point2);
        distanceTotal(i) = distance1(i) + distance2(i);
    end
    distanceMin = min(distanceTotal);
    distanceMinIndex = find(distanceTotal==distanceMin);
    ptIndex = k(distanceMinIndex);
    closestPoint = pointSet(ptIndex,:); % the point which is the closest to point 1 and 2
else % 3 points
    % Find the point which has the minimum sum of distances to point 1 ~ 3.
    ptsNum = size(pointSet,1); % unlabelled points number in the selected region
    centroidsTemp = pointSet; % temporary matrix used in the following processing
    k1 = []; k2 = []; k3 = [];
    for i = 1:6 % find 6 closest points to point 1 and 2 respectively
        if ptsNum == 0 % no point available, so stop the search
            break;
        end
        k1(i) = dsearchn(centroidsTemp,point1);
        k2(i) = dsearchn(centroidsTemp,point2);
        k3(i) = dsearchn(centroidsTemp,point3);
        centroidsTemp(k1(i),:) = -100; centroidsTemp(k2(i),:) = -100;
        centroidsTemp(k3(i),:) = -100; % set selected points to be invalid
        if k1(i) == k2(i) && k2(i) == k3(i)
            ptsNum = ptsNum - 1;
        elseif k1(i) == k2(i) || k2(i) == k3(i)
            ptsNum = ptsNum - 2;
        else
            ptsNum = ptsNum - 3;
        end
    end
    k0 = union(k1,k2); k = union(k0,k3); % indices of points to be judged
    ptsNumTBD = length(k);
    for i = 1:ptsNumTBD % calculate sums of distances
        distance1(i) = norm(pointSet(k(i),:)-point1);
        distance2(i) = norm(pointSet(k(i),:)-point2);
        distance3(i) = norm(pointSet(k(i),:)-point3);
        distanceTotal(i) = distance1(i) + distance2(i) + distance3(i);
    end
    distanceMin = min(distanceTotal);
    distanceMinIndex = find(distanceTotal==distanceMin);
    ptIndex = k(distanceMinIndex);
    closestPoint = pointSet(ptIndex,:); % the point which is the closest to point 1 and 2
end

end
