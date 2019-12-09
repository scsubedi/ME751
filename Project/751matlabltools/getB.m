function B = getB(p,abar)
% function that creates the matrix B(p,abar)
e = p(2:end);
B = 2*[(p(1)*eye(3)+tilda(e))*abar,  e*abar.' - (p(1)*eye(3)+tilda(e))*tilda(abar)];