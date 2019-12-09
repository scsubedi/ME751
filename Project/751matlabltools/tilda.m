function atilda = tilda(a)
% function that gives as output the cross-product matrix atilda for a
% vector a

atilda = [0,    -a(3),  a(2);
          a(3), 0,  -a(1);
          -a(2),    a(1),   0];