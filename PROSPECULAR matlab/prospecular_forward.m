% _______________________________________________________________________
% prospecular.m (August, 1 2023)
%
% This script used to forward simulate BRF, specular component and diffuse
% component.
% _______________________________________________________________________
% ***********************************************************************
% Xiao Li, Zhongqiu Suna, Shan Lua, Kenji Omasa,2023 
% PROSPECULAR: A model for simulating multi-angular spectral properties of 
% leaves by coupling PROSPECT with a specular function. Remote Sensing of Environment, 297, 113754
% ***********************************************************************
% Leaf BRF are calculated from 400 nm to 2500 nm (1 nm step) 

function [BRF,specular,diff]=prospecular_forward(a,wave,xr)
% a: biochemical parameters of PROSPECT model;
%       - n   = coefficient of leaf surface refractive index 
%       - rough   = roughness of leaf surface
%       - N   = leaf structure parameter
%       - Cab = chlorophyll a+b content in μg/cm2
%       - Car = carotenoids content in  μg/cm2
%       - Anth = Anthocyanin content in nmol/cm2
%       - Cbrown= brown pigments content in arbitrary units
%       - Cw  = equivalent water thickness in g/cm2 or cm
%       - Cm  = dry matter content in g/cm2
%       - Prot= protein content g/cm2
%       - CBC = Carbone-based constituents content in g/cm2 (cellulose, lignin, sugars...)
%
% xr: source zenith angle, viewing zenith angle and viewing azimuth angle;
% wave: wavelenth
%%  Diffuse component calculation
N=a(3);
Cab=a(4);
Ccx=a(5);
Can=a(6);
Cbp=a(7);
Cw=a(8);
Cdm=a(9);
Cpro=a(10);
Ccbc=a(11);
RT=prospect_PROdiff(N,Cab,Ccx,Can,Cbp,Cw,Cdm,Cpro,Ccbc);
R=RT(wave,2); %diffuse component
%%  specular component caliculation
% refractive index of leaf surface of 2009 Stucken
data    = data_external;
waven = data(wave,2);
rc(1,1:size(waven,1))=a(1).*waven';%refractive index
rc(2,:)=a(2).*ones(1,size(waven,1));%roughness

% illumination-viewing geometry
i=xr(1,:)';% source zenith angle, radian
e=xr(2,:)';% viewing zenith angle, radian
fai=xr(3,:)';% viewing azimuth angle, radian

% other required parameters
% Bousquet, L., Lachérade, S., Jacquemoud, S., & Moya, I., 2005.
% Leaf BRDF measurements and model for specular and diffuse components differentiation. Remote Sensing of Environment, 98, 201-211Leaf BRDF measurements and model for specular and diffuse components differentiation
ph=(acos(cos(i).*cos(e)+sin(i).*sin(e).*cos(fai)));
pha=(acos((cos(i)+cos(e))./(2.*cos((acos(cos(i).*cos(e)+sin(i).*sin(e).*cos(fai))/2)))));
E1=((2.*cos(pha).*cos(e))./cos(ph/2));
E2=((2.*cos(pha).*cos(i))./cos(ph/2));
G=(min(1,min(E1,E2)));
xn=[i,e,ph,pha,G];
g=sqrt(abs(rc(1,:).^2+cos(0.5.*xn(:,3)).^2-1));%g
fnq=0.5.*((g-cos(0.5.*xn(:,3)))./(g+cos(0.5.*xn(:,3)))).^2;%The first half of Fresnel factor
fnh=1+((cos(0.5.*xn(:,3)).*(g+cos(0.5.*xn(:,3)))-1)./(cos(0.5.*xn(:,3)).*(g-cos(0.5.*xn(:,3)))+1)).^2;%The second half of Fresnel factor
Fn=fnq.*fnh;%Fresnel factor
fzq=Fn.*xn(:,5);
fmq=4.*pi.*cos(xn(:,1)).*cos(xn(:,2));
fzh=exp(-tan(xn(:,4)).^2./rc(2,:).^2);
fmh=rc(2,:).*rc(2,:).*cos(xn(:,4)).^4;
brdfj=(fzq./fmq).*(fzh./fmh);%BRDF of specular component

%% output
BRF=R+brdfj'.*pi;%BRF=diffuse component +specular component
specular=brdfj'.*pi;%specular component
diff=R;%diffuse component
end
