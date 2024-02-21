function [n, m, q] = Noll2NMQ(nollIndex)

% Convert Noll index(j) to radial(n), angular(m) and orientation(q) indices
% for Standard Zernike Polynomials. Associations for q indices: 
%       q = +1 is cosine(angular term)
%       q = -1 is sine(angular term)
% note: code is vectorized over Noll indices, so the input can be an array.
% Revision History
% 2015-03-09 - Greg Smith - vectorized version of FindNMQ
% 2023-01-22 - Yiyang Huang - add notes

% Input validation.
if any(nollIndex <= 0)
    error('Noll2NMQ:invalidIndex', ...
        'Noll indices are required to be positive valued.');
end

% Compute n values.
n = ceil((sqrt(8*nollIndex+1)-3)/2); % from the relationship: [1+(n+1)]/2*(n+1) = j_max
    
% Compute m values.
k = nollIndex - n.*(n+1)./2; % the order of j in the sequence within level n
testOdd = find(mod(n, 2));
m = 2*floor(k/2); % even case of n
m(testOdd) = 2*floor((k(testOdd)-1)/2)+1; % odd case of n
    
% Compute q values (for discriminating even, odd and m==0 terms of j).
q = ones(size(n)); % even case of j, also indicates the cos term
q(m==0) = 0; % m==0 case
q(mod(nollIndex, 2) == 1 & m ~= 0) = -1; % odd case of j, also indicates the sin term

end
