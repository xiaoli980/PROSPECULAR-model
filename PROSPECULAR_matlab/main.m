% _______________________________________________________________________
% main.m
% version 1 (August, 1 2023)
% subroutines required: prospecular.m, prospect_PROdiff.m,
%                      calctav.m, dataSpec_PRO.m, data_external.m
%                      prospecular_forward.m
% _______________________________________________________________________
 
% This script allows to simulate multi-angular bidirectional reflectance factor (BRF) and 
% tp retrieves both leaf biochemical and surface structural parameters from BRF.
%
% Diffuse component (400-2500 nm) with 1 nm step is calculated using PROSPECT-PRO (Féret, Jean-Baptiste et al., 2021).
% The matlab codes of PROSPECT-PRO are available on the githab website of 
% Jean-Baptiste Feret at  https://gitlab.com/jbferet/prospect_pro_matlab.
%
% Specular component is calculated based on specular reflection function (Bousquet et al. 2005)
% _______________________________________________________________________
% 
% Xiao Li, Zhongqiu Sun, Shan Lu, Kenji Omasa, 
% PROSPECULAR: A model for simulating multi-angular spectral properties of 
% leaves by coupling PROSPECT with a specular function
%
% Féret, J.-B., Berger, K., de Boissieu, F., & Malenovský, Z., 2021. 
% PROSPECT-PRO for estimating content of nitrogen-containing leaf proteins 
% and other carbon-based constituents. Remote Sensing of Environment, 252, 112173
%
%Bousquet, L., Lachérade, S., Jacquemoud, S., & Moya, I., 2005.
% Leaf BRDF measurements and model for specular and diffuse components differentiation. Remote Sensing of Environment, 98, 201-211Leaf BRDF measurements and model for specular and diffuse components differentiation
% _______________________________________________________________________

clear 
clc
%% Spectral simulation

% load  structure and  biochemistry parameters
load('leaf_parameter.txt');%leaf parameter 
% n     = leaf_parameter(1);% coefficient of leaf surface refractive index
% rough    = leaf_parameter(2);% roughness of leaf surface
% N     = leaf_parameter(3);% leaf structure parameter
% Cab   = leaf_parameter(4);% chlorophyll a+b content in μg/cm2
% Car   = leaf_parameter(5);% carotenoids content in  μg/cm2
% Anth  = leaf_parameter(6); % Anthocyanin content in nmol/cm2
% Cbrown= leaf_parameter(7);% brown pigments content in arbitrary units
% Cw    = leaf_parameter(8);% equivalent water thickness in g/cm2 or cm
% Cm    = leaf_parameter(9);% dry matter content in g/cm2
% Prot  = leaf_parameter(10);protein content g/cm2
% CBC   = leaf_parameter(11); Carbone-based constituents content in g/cm2 (cellulose, lignin, sugars...)

% load illumination-viewing geometry
geo=load('geometry.txt');%leaf parameter
% SZA    = geo(1,:);% Source zenith angle,degree
% VZA    = geo(2,:);% Viewing zenith angle,degree
% VAA    = geo(3,:);% Viewing azimuth angle,degree
geometry=deg2rad(geo);% transfer degree into radian

% identification of wavelength
waveo=400:2500;
wave=waveo-399;% from 400 nm

% simulate
BRF=prospecular(leaf_parameter,wave,geometry);%simulate BRF

out=[geo;BRF]; % Correspond geometric parameters to BRF
xlswrite('leaf_spectrum.xlsx',out,'BRF')% output BRF

%% Model inversion

% identify parameter boundary and initial value
% P0=[n rough N  Cab  Car  Anth  Cbrown  Cw  Cm  Prot  CBC]
P0=[1.5  0.3  1.5	40	10	0.1	0	0.01	0	0.001	0.009];% initial value
LB=[1.0001  0.01  1	 0	0	0	0	0	0	0.0001	0.001]; %lower boundary
UB=[5	1  3.5	120	30	40	1	0.06	0	0.003	0.035];%upper boundary

options = optimset('Algorithm','trust-region-reflective','TolFun',1e-10);
sol= lsqcurvefit(@prospecular,P0,wave,BRF,LB,UB,options,geometry);% invert based on BRF,geometry parameters input as prior information

% reconstruct spectral 
[mnBRF,specular,diff]=prospecular_forward(sol,wave,geometry);% simulate specular and diffuse component with inverted parameters

% output
xlswrite('leaf_spectrum.xlsx',mnBRF,'Simulated_BRF')% output simulated BRF
xlswrite('leaf_spectrum.xlsx',specular,'specular')% output specular component
xlswrite('leaf_spectrum.xlsx',diff,'diff')% output diffuse component
%% Optinal :output data used for polar plot
data=BRF;
source=[geo(1,1);0];% SZA &SAA
wave=400:2500;%wavelength
plot=[[source;waveo'] [geo(2:3,:);data]];%Merging BRF, wavelengths and geometric parameters
p=mfilename('fullpath');%Find the current matlab file path (the path contains the matlab file name)
[filepath,~,~]=fileparts(p);%get path filepath
save([filepath,'\polar_plot\plot.mat'],'plot')  % save plot parameter to the 'polar plot' file 