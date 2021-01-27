function [pa, pdA, pdB] = prinangle(A, B)
%
% [pa, pdA, pdB] = prinangle(A, B)
%
% Principal angles and directions between column spaces of A and B
%
% INPUTS:
% A - defines p vectors in n-dimensional space (n x p)
% B - defines q vectors in n-dimensional space (n x q)
%
% OUTPUTS:
% pa  - principal angles 
% pdA - columns are principal directions in CS(A)
% pdB - columns are principal directions in CS(B)
%
% July 8, 2011 -- fixed bug for p=1 or q=1 case by adding 'econ' to svd 
% June 7, 2011 -- first created 
% @ 2011 Byron Yu -- byronyu@cmu.edu

if nargout == 1
  % More efficient if only principal angles are needed
  sv = svd(orth(A)' * orth(B));
  pa = acos(sv);
  
else
  Aorth = orth(A);
  Borth = orth(B);
  
  [UU, DD, VV] = svd(Aorth' * Borth, 'econ');
  
  pa  = acos(diag(DD));
  pdA = Aorth * UU;
  pdB = Borth * VV;
end