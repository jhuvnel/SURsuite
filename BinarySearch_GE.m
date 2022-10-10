% Return the index of the value in the sorted array ARY that is greater
% than or equal to VAL. This algorithm taken from Wikipedia!
%
%   https://en.wikipedia.org/wiki/Binary_search_algorithm
%
function L = BinarySearch_GE(ARY, VAL)
L = 1;
R = length(ARY)+1;
while L < R
   M = floor((L+R)/2);
%   fprintf('L %d, M %d, R %d\n', [L M R])
   if ARY(M) < VAL
      L = M+1;
   else 
      R = M;
   end
end

if L>length(ARY); L=[]; end

%fprintf('L %d, M %d, R %d\n', [L M R])
