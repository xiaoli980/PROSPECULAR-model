% plot_BRDF.m
% 
% Author: Laurent Bousquet
% Affiliation: Institut de Physique du Globe de Paris - University Paris Diderot
% Created: December 02, 2006
% Last modified: September 3, 2013
%
% This program plots the BRDF/BTDF of a plant leaf in a polar coordinate system
% It calls the function "FunPlotDir_pub"
clear
clc
load mnplot.mat %load plot data obtained from model forward
% filename=input('File name: ','s');
% mat2Bplot=load([filename,'.txt']);
mat2Bplot=plot;
[n,p]=size(mat2Bplot);
% nb_fig=input('Figue number: ');
 nb_fig=1;%figure 1
% lambda=input('Wavelength (400 nm < lambda < 2500 nm): '); % wavelength (nm)
lambda=670; %wavelength of plot 
line=lambda-mat2Bplot(3,1)+3;
theta_i=mat2Bplot(1,1); % illumination zenith angle (degrees)
phi_i=mat2Bplot(2,1);  % illumination azimuth angle (degrees)
theta_v=transpose(mat2Bplot(1,2:p)); % viewing zenith angles (degrees)
phi_v=transpose(mat2Bplot(2,2:p)); % viewing azimuth angles (degrees)
BRDF=transpose(mat2Bplot(line,2:p)); % leaf BRDF/BTDF at the selected wavelength
Vmin=0; % lower value of the scale bar (set to min of data if Vmin == -1)
Vmax=0; % upper value of the scale bar (set to max of data if Vmax == 0)

% FunPlotDir_pub_FunPlotDir_pub_10VZA_30VAA(nb_fig,lambda,theta_i,phi_i,theta_v,phi_v,BRDF,Vmin,Vmax);%original code
FunPlotDir_pub_10vza_30vaa(nb_fig,lambda,theta_i,phi_i,theta_v,phi_v,BRDF,Vmin,Vmax); % The viewing zenith angle (VAA) and viewing azimuth angle (VZA) intervals have been changed to 10° and 30°.

print(gcf,'-r600','-djpeg','mnBRF.jpeg');%save