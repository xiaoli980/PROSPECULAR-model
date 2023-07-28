% _______________________________________________________________________
% prospect_PROdiff.m
% This script allows to simulate diffuse component based on PROSPECT-PRO.
%
% prospect_PRO.m
% version 7.0 (January, 7th 2020)
% subroutines required: calctav.m, dataSpec_PRO.m
% _______________________________________________________________________
%
% Plant leaf reflectance and transmittance are calculated from 400 nm to
% 2500 nm (1 nm step) with the following parameters:
%
%       - N   = leaf structure parameter
%       - Cab = chlorophyll a+b content in µg/cm?
%       - Car = carotenoids content in µg/cm?
%       - Anth = Anthocyanin content in nmol/cm?
%       - Cbrown= brown pigments content in arbitrary units
%       - Cw  = equivalent water thickness in g/cm? or cm
%       - Cm  = dry matter content in g/cm?
%       - Prot= protein content g/cm?
%       - CBC = Carbone-based constituents content in g/cm? (cellulose, lignin, sugars...)
%
% Here are some examples observed during the LOPEX'93 experiment on
% fresh (F) and dry (D) leaves :
%
% ---------------------------------------------
%                N     Cab     Cw        Cm    
% ---------------------------------------------
% min          1.000    0.0  0.004000  0.001900
% max          3.000  100.0  0.040000  0.016500
% corn (F)     1.518   58.0  0.013100  0.003662
% rice (F)     2.275   23.7  0.007500  0.005811
% clover (F)   1.875   46.7  0.010000  0.003014
% laurel (F)   2.660   74.1  0.019900  0.013520
% ---------------------------------------------
% min          1.500    0.0  0.000063  0.0019
% max          3.600  100.0  0.000900  0.0165
% bamboo (D)   2.698   70.8  0.000117  0.009327
% lettuce (D)  2.107   35.2  0.000244  0.002250
% walnut (D)   2.656   62.8  0.000263  0.006573
% chestnut (D) 1.826   47.7  0.000307  0.004305
% ---------------------------------------------
% _______________________________________________________________________

% this code includes numerical optimizations proosed in the FLUSPECT code
% Authors: Wout Verhoef, Christiaan van der Tol (tol@itc.nl), Joris Timmermans, 
% Date: 2007
% Update from PROSPECT to FLUSPECT: January 2011 (CvdT)

function LRT=prospect_PROdiff(N,Cab,Car,Ant,Brown,Cw,Cm,Prot,CBC)
data    = dataSpec_PRO;
lambda  = data(:,1);    nr      = data(:,2);
Kab     = data(:,3);    Kcar    = data(:,4);
Kant    = data(:,5);    Kbrown  = data(:,6);
Kw      = data(:,7);    Km      = data(:,8);
Kprot   = data(:,9);    Kcbc    = data(:,10);

Kall=(Cab*Kab+Car*Kcar+Ant*Kant+Brown*Kbrown+Cw*Kw+Cm*Km+Prot*Kprot+CBC*Kcbc)/N;

j           = find(Kall>0);               % Non-conservative scattering (normal case)
t1          = (1-Kall).*exp(-Kall);
t2          = Kall.^2.*expint(Kall);
tau         = ones(size(t1));
tau(j)      = t1(j)+t2(j);

% ***********************************************************************
% reflectance and transmittance of one layer
% ***********************************************************************
% Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969),
% Interaction of isotropic ligth with a compact plant leaf, J. Opt.
% Soc. Am., 59(10):1376-1379.
% ***********************************************************************
% reflectivity and transmissivity at the interface
%-------------------------------------------------
talf        = calctav(40,nr);
ralf        = 1-talf;
t12         = calctav(90,nr);
r12         = 1-t12;
t21         = t12./(nr.^2);
r21         = 1-t21;

% top surface side
denom       = 1-r21.*r21.*tau.^2;
Ta          = talf.*tau.*t21./denom;
Ra          = ralf+r21.*tau.*Ta;

% bottom surface side
t           = t12.*tau.*t21./denom;
r           = r12+r21.*tau.*t;

% ***********************************************************************
% reflectance and transmittance of N layers
% Stokes equations to compute properties of next N-1 layers (N real)
% Normal case
% ***********************************************************************
% Stokes G.G. (1862), On the intensity of the light reflected from
% or transmitted through a pile of plates, Proc. Roy. Soc. Lond.,
% 11:545-556.
% ***********************************************************************

D           = sqrt((1+r+t).*(1+r-t).*(1-r+t).*(1-r-t));
rq          = r.^2;
tq          = t.^2;
a           = (1+rq-tq+D)./(2*r);
b           = (1-rq+tq+D)./(2*t);

bNm1        = b.^(N-1);                  %
bN2         = bNm1.^2;
a2          = a.^2;
denom       = a2.*bN2-1;
Rsub        = a.*(bN2-1)./denom;
Tsub        = bNm1.*(a2-1)./denom;

%			Case of zero absorption
j           = find(r+t >= 1);
Tsub(j)     = t(j)./(t(j)+(1-t(j))*(N-1));
Rsub(j)	    = 1-Tsub(j);

% Reflectance of the leaf: combine top layer with next N-1 layers
denom       = 1-Rsub.*r;
refl        = Ra-ralf+Ta.*Rsub.*t./denom;
LRT = [lambda refl];

