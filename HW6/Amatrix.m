
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtains the A matrix (Direction cosine matrix or Rotation transformation
% matrix from euler parameters). 
% Input: Euler Parameters (e0 (1 x 1), e(3 x 1)
% Output: A(3 x 3) matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = Amatrix(e0, e)

G = [-e -tilda(e)+e0*eye(3)];
E = [-e, tilda(e) + e0*eye(3)];
A = E*G';

end
