function zernikeVal = Zernike(rho, theta, zernikeTerms, options)


% ZERNIKE computes Zernike values at specific position and for specific terms
% INPUT:
%   rho: radial position of Zernike coefficient (normalized to unit radius)
%   theta: angular position of Zernike coefficient [rad]
%   zernikeTerms: 1 or more Zernike terms to calculate (Noll index)
%   options: structure with the following optional elements
%       .largeRho: If set to 'keep', then radial positions outside unit 
%           circle will be computed. If set to 'ignore', radial positions 
%           outside unit circle will have Zernike values set to NaN. [default = 'ignore']
% OUTPUT:
%   zernikeVal: array of Zernike values. Each rho and theta will have an 
%       array of Zernike values. For instance, if rho and theta are 2D
%       arrays, and ZernikeTerms is 1D array of numbers, then zernikeVal
%       will be a 3D array of dimensions [size(rho,1), size(rho,2),
%       numel(zernikeTerms)]. Values are normalized to unit rms.
% DEPENDENCIES:
% - Noll2NMQ
% NOTES:
% - rho and theta can be single number, 1D array, or 2D array of values but both must have the same dimensions
% - zernikeTerms can be a single number or 1D array
% - beyond Noll index of approx 1,000 the results may have calculation artifacts on the order of the signal strength if rho is near 1
% - beyond Noll index of approx 14,700 the result will contain invalid values because intermediate calculations exceed machine storage limitations
% Revision History:
% 2015-03-09 - Greg Smith - adapted from ZStdPA and related code
% 2017-12-18 - Greg Smith - add optional parameter to allow computations outside unit circle
% 2023-01-23 - Yiyang Huang - add notes


% Sanitize input.
if nargin < 4 || ~isfield(options,'largeRho'), options.largeRho = 'ignore'; end
if size(rho) ~= size(theta)
    error('Zernike:sizeMismatch','Rho and theta must have the same dimensions.');
end
% How to handle large radial positions outside the unit circle.
switch lower(options.largeRho)
    case 'ignore' % ignore points outside the unit circle
        nanFlag = true; 
    case 'keep' % still compute points outside the unit circle
        nanFlag = false;
    otherwise
        error('Zernike:badOption',['''largeRho'' optional input must ' ...
            'be set to ''ignore'' or ''keep''']);
end


% Convert from Noll indices to subscript/superscript notations.
[n, m, q] = Noll2NMQ(zernikeTerms);
% Radial dependence.
radial = Radial(n, m, rho, nanFlag); % 3D radial dependence (the 3rd dimension shows the result of every polynomial)
% Angular dependence.
angleFactor = ApplyRepeated(theta, m, ones(size(m))); % 3D angular data
angular = ones([size(theta),numel(q)]); % default
angular(:,:,q==1) = sqrt(2)*cos(angleFactor(:,:,q==1)); % q == +1
angular(:,:,q==-1) = sqrt(2)*sin(angleFactor(:,:,q==-1)); % q == -1
% Combine to output.
scaleFactor = permute(repmat(sqrt(n(:)+1),[1,size(theta)]),[2,3,1]);
zernikeVal = scaleFactor.*radial.*angular; % every layer of the 3D data represents a single Zernike polynomial


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radial = Radial(n, m, rho, nanFlag)
% Radial dependence of the Zernike standard polynomials.

% Initialize output.
radial = zeros([size(rho),numel(n)]);
% Compute the radial term of Zernike polynomials.
sign = -1;
sMax = (n-m)/2; % maximum iteration count for all polynomials
for s = 0:max(sMax) % summation process - every layer works for all polynomials
    sign = -1 * sign; % alternate sign for each iteration of summation
    term = FactorialTerm(n, m, s*ones(size(n))); % 1D factorial terms for all polynomials on a single layer
    addition = ApplyRepeated(rho, sign*term, (n-2*s)); % 3D radial dependence on a single layer (the 3rd dimension size equals to the polynomials length)
    radial(:,:,s<=sMax) = radial(:,:,s<=sMax) + addition(:,:,s<=sMax);
end
% Set any values with rho > 1 (outside the unit circle) to nan.
if nanFlag == true
    rhoCube = repmat(rho,[1,1,numel(n)]);
    radial(rhoCube > 1) = nan;
end

end


function value = FactorialTerm(n, m, s)
% Compute factorial term of Zernike radial function.
% n, m, and s must have same dimensions (vectorized over n and m).
% note: if n < s or (n-m)/2-s < 0 then result is nan.

% Check if any term will exceed factorial limit. The limit from MATLAB 
% helps assuming double precision numbers (no quad precision available).
factorialLimit = 170;
factor1 = n-s;
factor2 = s;
factor3 = (n+m)/2-s;
factor4 = (n-m)/2-s;
if any([factor1, factor2, factor3, factor4] > factorialLimit)
    warning('Zernike:radial', ...
        ['This calculation contains factorial terms which exceed the ' ...
        'limits of machine number storage. Results may be invalid.']);
end
    
% Make note of negative factors and set to zero for factorial calculation.
neg = factor1<0 | factor2<0 | factor3<0 | factor4<0;
factor1(neg) = 0; factor2(neg) = 0; factor3(neg) = 0; factor4(neg) = 0;
% Compute value and return after setting negative factorial terms to nan.
value = factorial(factor1)./ factorial(factor2)./ factorial(factor3)./ ...
    factorial(factor4);
value(neg) = nan; % usually means current s doesn't work for n and m

end


function cube = ApplyRepeated(matrix, scale, pow)
% Given 1D array of scale factors and powers, compute 3D array of matrices
% scale and power must be 1D arrays with the same number of elements.
% note: this is specifically helpful for radial terms.

% Change dimensionality of scale and power inputs to be like 3D output.
scale = reshape(scale,1,1,numel(scale));
pow = reshape(pow,1,1,numel(pow));
    
% Define function and apply calculation. This is slowest part. Further 
% optimization must occur at higher algorithm level.
funcHandle = @(x,y)x.*matrix.^y;
cube = bsxfun(funcHandle,scale,pow);

end
