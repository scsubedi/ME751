function varargout = getd(gcon_inp,varargin)
% function that calculates phi, nu,gamma, dphi_dr, dphi_dp for GCon D


%%%%%%%%%%%% enter vararg conditions here %%%%%%%%%%%%%%
% extracting input variables from structure
s1_p = gcon_inp.s1_p;
p1 = gcon_inp.p1;
p1dot = gcon_inp.p1dot;
s2_q = gcon_inp.s2_q;
p2 = gcon_inp.p2;
p2dot = gcon_inp.p2dot;
r1 = gcon_inp.r1;
r1dot = gcon_inp.r1dot;
r2 = gcon_inp.r2;
r2dot = gcon_inp.r2dot;
ft = gcon_inp.ft;
dft = gcon_inp.dft;
ddft = gcon_inp.ddft;
i = gcon_inp.i;
j = gcon_inp.j;

A1 = getA(p1);
A2 = getA(p2);

% calculating all the B matrices required
B_p2s2q = getB(p2,s2_q);
B_p2dots2q = getB(p2dot,s2_q);   
B_p1s1p = getB(p1,s1_p);
B_p1dots1p = getB(p1dot,s1_p);

d12 = (r2 + A2*s2_q - r1 - A1*s1_p);
d12dot = r2dot + B_p2s2q*p2dot - r1dot - B_p1s1p*p1dot;
phi = d12.'*d12 - ft;

nu = dft;

gamma = -2*d12.'*B_p2dots2q*p2dot + 2*d12.'*B_p1dots1p*p1dot - 2*(d12dot.'*d12dot) + ddft;

ddr1 = -2*d12.'*eye(3);               % DP1 has no r terms
ddr2 = 2*d12.'*eye(3);
ddp1 = -2*d12.'*B_p1s1p;
ddp2 = 2*d12.'*B_p2s2q;

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