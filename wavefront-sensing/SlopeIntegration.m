function surface = SlopeIntegration(slopeX, slopeY, type)


% Integrate surface gradient using matrix solution given in Southwell paper.
% The sampled slope data are assumed to align with the surface data. The
% calculation leverages Matlab's ability to handle large sparse matrices to
% directly solve the linear system of equations. The calculations closely
% follow Southwell's paper. See reference for details.
% INPUT:
%   slopeX: 2D array of numbers representing the horizontal gradient of the
%       surface. NaN values denote missing data.
%   slopeY: 2D array of numbers representing the vertical gradient of the
%       surface. NaN values denote missing data. Must have same dimensions 
%       as slopeX.
%   type: string indicator for the type of integration kernel to use. The 
%       following options are available:  [default = 'Southwell']
%       type = 'GradientInverse'
%           This type uses a kernel which matches the MATLAB 'gradient'
%           function. The gradient is computed as the symmetric
%           difference of the +1 pt and the -1 pt. For slope data
%           without NaN values, it will exactly undo the gradient
%           function. However because this function skips points, it
%           can produce checkerboard artifacts in the presence of sharp
%           edges.
%       type = 'Southwell'
%           This type uses the kernel given by equation 5 in the
%           reference paper. The gradient is the difference between
%           current point and the +1 pt. Because the kernel is
%           different from the MATLAB 'gradient' function, it cannot
%           undo that calculation. However the kernel does not skip
%           points, and therefore produces smoother output. [default]
% OUTPUT:
%   surface: 2D array of numbers representing the integrated slope data.
%   For every non-NaN gradient point, there will be a surface point.
% REFERENCES:
%   Southwell, W. "Wave-front estimation from wave-front slope 
%       measurements," J. Opt. Soc. Am.  70, 998-1006 (1980). 
%       https://doi.org/10.1364/JOSA.70.000998
% NOTES:
%   - Calculation assumes unit distance between points. One way to handle 
%       arbitrary distance is to multiply slope data by gradient of
%       position data.
%   - There is an arbitrary constant of integration that is undefined. This 
%       code sets the mean of all non-NaN values to zero prior to returning.
% HISTORY:
% 2016/05/21 - Isaac Trumper
% 2017/01/11 - Greg Smith - removed SAGUARO wrapper and units. Add padding for edge pixels. Create subroutines.
% 2017/02/20 - Greg Smith - modify to be inverse of built-in 'gradient' function
% 2017/03/01 - Greg Smith - add ability to switch between Gradient Inverse and Southwell kernel
% 2017/10/06 - Greg Smith - move 'type' into options structure, and add option for filling data
% 2017/10/10 - Greg Smith - Eliminate data filling option. Eliminate NaNs from sparse matrix instead.
% 2017/10/24 - Greg Smith - remove NaN data to reduce memory usage and be more explicit about removing bad rows from matrix
% 2021/04/26 - Stephanie Meyen - included threshold for bad data to define
% any data <threshold as bad data as well, not just NaNs.% edge problems on GMT3, needs to be revised every measurement cycle. lines
% 83-90, 105-107, 116-122


    % default inputs
    if nargin < 3, type = 'Southwell'; end

    % validate input (options are validated where used)
    if numel(size(slopeX)) ~= 2 || numel(size(slopeY)) ~= 2
        error('SlopeIntegration:wrongDim',...
            'SlopeX and slopeY must both be 2 dimensional arrays. Input slopeX has %i dimensions and slopeY has %i.',...
            numel(size(slopeX)),numel(size(slopeY))...
        );
    end
    if any(size(slopeX) ~= size(slopeY))
        error('SlopeIntegration:badDim',...
            'Dimensions of slopeX [%i, %i] must be same as slopeY [%i, %i].',...
            size(slopeX,1),size(slopeX,2),size(slopeY,1),size(slopeY,2)...
        );
    end


    % Create data mask for output
    % If either input is missing a point, it should be considered unreliable.
   %  badSurface = isnan(slopeX) | isnan(slopeY); % as it was before SM
    % change
    thresholdOn = 0; %use 0 for normal processing. Use 1 for processing that throws out data < threshold
    threshold = -5*10^-5; %added by SM (-3*10^-5 for GMT3)
    
    if thresholdOn == 1
    badSurface = isnan(slopeX) | isnan(slopeY)| slopeX(:,:)< threshold | slopeY(:,:)< threshold ; %added by SM
    else
    badSurface = isnan(slopeX) | isnan(slopeY); % as it was before SM    
    end
    
    % type-specific pre-processing
    switch lower(type)
        case 'southwell'
            % Southwell integration needs to be shifted by half-pixel to have surface pixels match slope pixels.
            [slopeX, slopeY] = SouthwellShift(slopeX, slopeY);
    end


    % Concatenate slope data and remove NaNs to minimize memory usage
    % Calculation will ignore missing data, so the actual value does not affect the output.
    data = [slopeX(:);slopeY(:)];
    data(isnan(data)) = 0; 
    
    if thresholdOn == 1
    data(data(:)< threshold) = 0; % added by SM 
    end
    
    % Size of the sparse matrix (M x N on each side)
    M = size(slopeX,1);
    N = size(slopeX,2);
    
     
    % Create sparse matrix that describes the conversion from surface to gradient.
  %  sparseMatrix = CreateMatrix(~isnan(slopeX), ~isnan(slopeY), type);
    if thresholdOn == 1
      goodSurfaceX = ~isnan(slopeX) & slopeX(:,:)> threshold; %added by SM
      goodSurfaceY = ~isnan(slopeY) & slopeY(:,:)> threshold; %added by SM    
    sparseMatrix = CreateMatrix(goodSurfaceX, goodSurfaceY, type);
    else 
    sparseMatrix = CreateMatrix(~isnan(slopeX), ~isnan(slopeY), type); %was before SM
    end
    
    % Solve the linear system of equations
    % During calculation, a warning is generated about a rank-deficient
    % matrix because the constant of integration is not defined. This is
    % discussed in Southwell's document (see reference). One can remedy
    % this by adding another condition to the matrix equation, but practice
    % shows that doing so increases time and memory by orders of magnitude.
    % It is easier to manually remove the mean later and temporarily turn
    % off rank-deficient warning during the calculation. Additionally, 
    % locations with missing data also lead to undefined surface values.
    warningStruct = warning;
    warning('off','MATLAB:rankDeficientMatrix');
        surface = sparseMatrix\data;
    warning(warningStruct);

    % Reshape the sag vector, which was computed in column major format, in
    % to a 2D array that represents the actual surface.
    surface = reshape(surface,M,N);

    % re-apply the mask
    surface(badSurface) = NaN;

    % Manually remove the mean value. This is equivalent to setting piston
    % to zero and is the arbitrary constant of integration.
    surface = surface - mean(surface(:),'omitnan');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sparseMatrix = CreateMatrix(goodDataX, goodDataY, type)
% This function uses pre-defined arrays to match the 'gradient' function output
% If the 'gradient' changes, these arrays will need to change.
%
% INPUT:
%   goodDataX, goodDataY: 2D array of 1s and 0s where 1s represent data to 
%       be processed and 0s represent bad data (NaN values). Array size
%       must match dimensions of gradient image and must be same for X & Y.
%   type: string indicator for the type of integration kernel to use. ('GradientInverse', or 'Southwell')
%
% OUTPUT:
%   sparse Matrix describing conversion from 2D surface array to two 2D 
%   gradient arrays. If slopeSize is [M,N], then output will be sparse 
%   array of size [2*M*N, M*N]. The number of rows is double because
%   Y-gradient comes after the X-gradient. The resultant matrix will 
%   convert a 1D flattened surface array [surface(:)] to a single list 
%   which can be represented as [gradX(:);gradY(:)];
%
% NOTE:
%   Empty rows in the sparse matrix will be associated with data which is 
%   not good. This makes it so the integrated results do not depend on 
%   those values.

% HISTORY:
% 2017-10-24 - Greg Smith - adjust to allow for data masking

    % Define convolution kernel. This is effectively stating how the 
    % equivalent gradient function handles adjacent pixels. The kernel is
    % defined by three points: [-1 neighbor, current pixel, +1 neighbor]
    
    
    switch lower(type)
        case 'gradientinverse'
            % For most cases, use symmetric separation (2 pixels apart), but for edges, use 1 pixel apart.

            % convolution kernel for interior points
            interiorX = [-0.5, 0, 0.5];
            interiorY = [-0.5; 0; 0.5];

            % convolution kernel for top points
            topX = [-0.5, 0, 0.5];
            topY = [0; -1; 1];

            % convolution kernel for right points
            rightX = [-1, 1, 0];
            rightY = [-0.5; 0; 0.5];

            % convolution kernel for bottom points
            bottomX = [-0.5, 0, 0.5];
            bottomY = [-1; 1; 0];

            % convolution kernel for left points
            leftX = [0, -1, 1];
            leftY = [-0.5; 0; 0.5];
        case 'southwell'
            % Use current pixel and +1 neighbor for calculation. (see right-side of reference equations 9-10)
            % Southwell doesn't specify how to handle edge points, so shift things over as an approximation

            % convolution kernel for interior points
            interiorX = [0, -1, 1];
            interiorY = [0; -1; 1];

            % convolution kernel for top points
            topX = [0, -1, 1];
            topY = [0; -1; 1];

            % convolution kernel for right points
            rightX = [-1, 1, 0];        % special
            rightY = [0; -1; 1];

            % convolution kernel for bottom points
            bottomX = [0, -1, 1];
            bottomY = [-1; 1; 0];       % special

            % convolution kernel for left points
            leftX = [0, -1, 1];
            leftY = [0; -1; 1];
        otherwise
            error('SlopeIntegration:badType',...
                '''options.type'' input (''%s'') not recognized by SlopeIntegration function. Allowed types are: ''GradientInverse'' and ''Southwell''.', type);
    end

    % get some array indices to be used for larger (sparse) matrix (creates
    % MxN matrix with index numbering from upper left corner starting at 1
    % and going down all lines and contiune at top column by column
    
    
    [M, N] = size(goodDataX);
    indexLin = reshape(1:M*N,M,N);


    % Determine sparse matrix indices and values for every point
        % interior points
        selectPts = indexLin(2:end-1,2:end-1);
        indexArrayX = CreateIndexArray(selectPts(:), [M,0], interiorX, goodDataX(2:end-1,2:end-1));
        indexArrayY = CreateIndexArray(selectPts(:), [0,1], interiorY, goodDataY(2:end-1,2:end-1));

        % top points (except corner)
        selectPts = indexLin(1,2:end-1);
        indexArrayX = [indexArrayX; CreateIndexArray(selectPts(:), [M,0], topX, goodDataX(1,2:end-1))];
        indexArrayY = [indexArrayY; CreateIndexArray(selectPts(:), [0,1], topY, goodDataY(1,2:end-1))];

        % right points (except corner)
        selectPts = indexLin(2:end-1,end);
        indexArrayX = [indexArrayX; CreateIndexArray(selectPts(:), [M,0], rightX, goodDataX(2:end-1,end))];
        indexArrayY = [indexArrayY; CreateIndexArray(selectPts(:), [0,1], rightY, goodDataY(2:end-1,end))];

        % bottom points (except corner)
        selectPts = indexLin(end,2:end-1);
        indexArrayX = [indexArrayX; CreateIndexArray(selectPts(:), [M,0], bottomX, goodDataX(end,2:end-1))];
        indexArrayY = [indexArrayY; CreateIndexArray(selectPts(:), [0,1], bottomY, goodDataY(end,2:end-1))];

        % left points (except corner)
        selectPts = indexLin(2:end-1,1);
        indexArrayX = [indexArrayX; CreateIndexArray(selectPts(:), [M,0], leftX, goodDataX(2:end-1,1))];
        indexArrayY = [indexArrayY; CreateIndexArray(selectPts(:), [0,1], leftY, goodDataY(2:end-1,1))];

        % manually determine corners since there are only 4 of them
        indexArrayX = [indexArrayX; ...
            [indexLin(1,1), indexLin(1,1), leftX(2).*goodDataX(1,1)];               % upper-left
            [indexLin(1,1), indexLin(1,1)+M, leftX(3).*goodDataX(1,1)];...          % upper-left
            [indexLin(M,1), indexLin(M,1), leftX(2).*goodDataX(end,1)];             % lower-left
            [indexLin(M,1), indexLin(M,1)+M, leftX(3).*goodDataX(end,1)];...        % lower-left
            [indexLin(M,N), indexLin(M,N), rightX(2).*goodDataX(1,end)];            % lower-right
            [indexLin(M,N), indexLin(M,N)-M, rightX(1).*goodDataX(1,end)];...       % lower-right
            [indexLin(1,N), indexLin(1,N), rightX(2).*goodDataX(end,end)];          % upper-right
            [indexLin(1,N), indexLin(1,N)-M, rightX(1).*goodDataX(end,end)]];       % upper-right
        indexArrayY = [indexArrayY; ...
            [indexLin(1,1), indexLin(1,1), topY(2).*goodDataY(1,1)];                % upper-left
            [indexLin(1,1), indexLin(1,1)+1, topY(3).*goodDataY(1,1)];...           % upper-left
            [indexLin(M,1), indexLin(M,1), bottomY(2).*goodDataY(end,1)];           % lower-left
            [indexLin(M,1), indexLin(M,1)-1, bottomY(1).*goodDataY(end,1)];...      % lower-left
            [indexLin(M,N), indexLin(M,N), bottomY(2).*goodDataY(1,end)];           % lower-right
            [indexLin(M,N), indexLin(M,N)-1, bottomY(1).*goodDataY(1,end)];...      % lower-right
            [indexLin(1,N), indexLin(1,N), topY(2).*goodDataY(end,end)];            % upper-right
            [indexLin(1,N), indexLin(1,N)+1, topY(3).*goodDataY(end,end)]];         % upper-right

    % Combine index arrays
    % This requires offsetting indexArrayY to be after indexArrayX
    indexArrayY(:,1) = indexArrayY(:,1) + M*N;
    indexArray = [indexArrayX;indexArrayY];

    % Last chance to throw away any data with zeros for the value.
    % Most likely, this only cleans up the manually-entered corner values.
     indexArray = indexArray(abs(indexArray(:,3)) > 10^-10,:);
 

    % Turn indices and values into a large, sparse matrix
    sparseMatrix = sparse(indexArray(:,1), indexArray(:,2), indexArray(:,3), 2*M*N, M*N);
end



function indexArray = CreateIndexArray(indices, shift, template, goodData)
% This is a way to automate the index handling for the sparse matrix
%
% INPUT:
%   indices: 1D array of matrix indices at center of pixel
%   shift: 2-element array specifying how much offset to handle +/-
%   template: 3-element array with the -shift, center, and +shift values
%   goodData: Array of 1s (good) and 0s (bad) describing where the data is 
%       useful for processing. Must have same number of elements as 
%       'indices' input.
%
% OUTPUT:
%   Nx3 array where columns are [row index, col index, value]. The length 
%   N is the number of non-zero elements in template times the number of 
%   indices.
%
% NOTE: template values of zero are not returned.

% HISTORY:
% 2017-10-24 - Greg Smith - added goodData input

    % handle points at specific index
    indexCenter = [indices, indices, repmat(template(2),size(indices)).*goodData(:)];

    % handle points shifted in X
    indexPlusX = [];
    indexMinusX = [];
    if shift(1) > 0
        indexPlusX = [indices, indices+shift(1), repmat(template(3),size(indices(:))).*goodData(:)];
        indexMinusX = [indices, indices-shift(1), repmat(template(1),size(indices(:))).*goodData(:)];
    end

    % handle points shifted in Y
    indexPlusY = [];
    indexMinusY = [];
    if shift(2) > 0
        indexPlusY = [indices, indices+shift(2), repmat(template(3),size(indices(:))).*goodData(:)];
        indexMinusY = [indices, indices-shift(2), repmat(template(1),size(indices(:))).*goodData(:)];
    end

    % combine for output and remove any points with zero value
    % (points with zero value don't do anything in sparse matrix)
    indexArray = [indexCenter; indexPlusX; indexMinusX; indexPlusY; indexMinusY];
    indexArray = indexArray(abs(indexArray(:,3)) > 10^-10,:);
    %indexArray = indexArray(abs(indexArray(:,3)) > threshold,:); %changed by SM
end



function [slopeXOut, slopeYOut] = SouthwellShift(slopeXIn, slopeYIn)
% Shift coordinates by half-pixel to make surface map match slope map.
% 
% The difference between the coincident surface/slope grid (see reference
% equation 9-10) and the Hudgin-aligned grid (see reference equation 5-6)
% is the same as the difference between figure 1a and figure 1b in the
% reference. The right-hand sides of the equations are the same, but the
% left-hand sides are different.
% 
% Assuming the integration carries out the right-hand side of equations
% 9-10, The only step is to re-interpret the slope data based on the left
% side of equations 9-10. The purpose of this function is to perform a
% half-pixel shifting so that after the half-pixel shift during integration,
% then the output will be same pixel as input.
%
% INPUT:
%   slopeXIn: 2D array representing the x-slope data to be integrated.
%   slopeYIn: 2D array representing the y-slope data to be integrated.
%
% OUTPUT:
%   slopeXOut: shifted slope such that after integration, the surface point 
%       will align with the original x-slope data.
%   slopeYOut: same as slopeX, but for the y-slope data.
%
% REFERENCES:
%   Southwell, W. "Wave-front estimation from wave-front slope 
%       measurements," J. Opt. Soc. Am.  70, 998-1006 (1980). 
%       https://doi.org/10.1364/JOSA.70.000998

% HISTORY:
% 2017-10-09 - Greg Smith - initial implementation.

    slopeXOut = slopeXIn;
    slopeXOut(:,1:end-1) = (slopeXIn(:,2:end) + slopeXIn(:,1:end-1))/2;
    slopeXOut(:,end) = slopeXOut(:,end-1); 

    slopeYOut = slopeYIn;
    slopeYOut(1:end-1,:) = (slopeYIn(2:end,:) + slopeYIn(1:end-1,:))/2;
    slopeYOut(end,:) = slopeYOut(end-1,:); 

end
