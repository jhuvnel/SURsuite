% Return the index of the value in the sorted array ARY that is LESS THAN
% or equal to VAL. This algorithm taken from Wikipedia!
%
%   https://en.wikipedia.org/wiki/Binary_search_algorithm
%
function L = BinarySearch_LE(ARY, VAL)
L = 1;
R = length(ARY)+1;
while L < R
   M = floor((L+R)/2);
%   fprintf('L %d, M %d, R %d\n', [L M R])
   if ARY(M) <= VAL  % Note for a LE search, this must be <= !!!
      L = M+1;
   else 
      R = M;
   end
end

L=L-1;
if L<1; L=[]; end

%fprintf('L %d, M %d, R %d\n', [L M R])
