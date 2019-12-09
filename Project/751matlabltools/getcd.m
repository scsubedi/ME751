function varargout = getcd(gcon_inp,varargin)
% function that calculates phi, nu,gamma, dphi_dr, dphi_dp for GCon CD


%%%%%%%%%%%% enter vararg conditions here %%%%%%%%%%%%%%
% extracting input variables from structure
s1_p = gcon_inp.s1_p;       % given vector in local coordinates, i.e. a1bar
p1 = gcon_inp.p1;
p1dot = gcon_inp.p1dot;
s2_q = gcon_inp.s2_q;
p2 = gcon_inp.p2;
p2dot = gcon_inp.p2dot;
c = gcon_inp.c;
r1 = gcon_inp.r1;
r2 = gcon_inp.r2;
ft = gcon_inp.ft;
dft = gcon_inp.dft;
ddft = gcon_inp.ddft;
i = gcon_inp.i;
j = gcon_inp.j;

A1 = getA(p1);
    s1p = A1*s1_p;        % finding vector in global coordinates
A2 = getA(p2);
    s2q = A2*s2_q;

% calculating all the B matrices required
B_p2 = getB(p2,s2_q);
B_p2dot = getB(p2dot,s2_q);   
B_p1 = getB(p1,s1_p);
B_p1dot = getB(p1dot,s1_p);
    
phi = c.'*(r2 + s2q - r1 - s1p) - ft;

nu = dft;

gamma = c.'*B_p1dot*p1dot - c.'*B_p2dot*p2dot + ddft;

ddr1 = -c.'*eye(3);               % DP1 has no r terms
ddr2 = c.'*eye(3);
ddp1 = -c.'*B_p1;
ddp2 = c.'*B_p2;

% checking if either bodies are ground
if i == 0
    dphi_dr = ddr2;
    dphi_dp = ddp2;
elseif j == 0
    dphi_dr = ddr1;
    dphi_dp = ddp1;
else
    dphi_dr = [ddr1,ddr2];
    dphi_dp = [ddp1,ddp2];
end

arg = 1;
if nargin > 1
    if any(strcmp('phi',varargin{1}))
        varargout{arg} = phi;
        arg = arg+1;
    end
    if any(strcmp('nu',varargin{1}))
        varargout{arg} = nu;
        arg = arg+1;
    end
    if any(strcmp('gamma',varargin{1}))
        varargout{arg} = gamma;
        arg = arg+1;
    end
    if any(strcmp('dphi_dr',varargin{1}))
        varargout{arg} = dphi_dr;
        arg = arg+1;
    end
    if any(strcmp('dphi_dp',varargin{1}))
        varargout{arg} = dphi_dp;
    end
else
    varargout = {phi,nu,gamma,dphi_dr,dphi_dp};
end