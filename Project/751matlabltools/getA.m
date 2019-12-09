function A = getA(p)
% function that gives as output A for a set of euler parameters p
e0 = p(1);
e1 = p(2);
e2 = p(3);
e3 = p(4);

A = 2*[e0^2+e1^2-0.5,   e1*e2-e0*e3,   e1*e3+e0*e2;
    e1*e2+e0*e3,    e0^2+e2^2-0.5,  e2*e3-e0*e1;
    e1*e3-e0*e2,    e2*e3+e0*e1,    e0^2+e3^2-0.5];