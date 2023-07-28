function FunPlotDir_pub_10vza_30vaa(nb_fig,lambda,thetai,phii,theta,phi,V,Vmin,Vmax)

% Compared to the original 'FunPlotDir_pub', the main differences of this version are as follows:
% 1. The VAA (Viewing Zenith Angle) and VZA (Viewing Azimuth Angle) intervals have been changed to 10бу and 30бу, respectively.
% 2. The viewing zenith angle labels have been added.
% 3. The sampling method has been changed from 'cubic' to 'natural'.
% 4. The position of the scale bar has been adjusted to avoid overlap with the polar plot.
%
% Modified by Xiao Li (August 1, 2013)
% Northeast normal university
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original information
% Authors: Laurent Bousquet and Cedric Bacour
% Affiliation: Institut de Physique du Globe de Paris - University Paris Diderot
% Created: November 22, 2006
% Last modified: September 3, 2013
%
% This function plots bidirectional data
% Depending on the option, it may call the function "geostat_gonio.m"

% input parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nb_fig = scalar, number of the current figure
% lambda = scalar, wavelength (nm)
% thetai = scalar, illumination zenith angle (degrees)
% phii = scalar, illumination azimuth angle (degrees)
% theta = n-by-1 vector, viewing zenith angles (degrees)
% phi = n-by-1 vector, viewing azimuth angles (degrees)
% V = vector N lines * 1 column, signal
% Vmin = scalar, lower value of the scale bar (set to min of data if Vmin == -1)
% Vmax = scalar, upper value of scale bar (set to max of data if Vmax == 0)

% test set
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 20 ;
% lambda = 590 ;
% thetai = 40 ;
% phii = 0 ;
% theta = [0:90/(N-1):90] ;
% phi = [0:360/(N-1):360] ;
% V = [0:1/(N-1):1] ;
% Vmin = 0 ;
% Vmax = 1 ;

% 1. Prepare the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2r = (pi/180) ;
%figure(nb_fig) ;
hold on ;

% 2. Plot features independent of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1. Window
h=plot([-90 -90 90 90],[-90 90 90 -90]);
axis image;
delete(h) ;
% 2.2. Outer circle
th=0:pi/50:2*pi;
xunit=90*cos(th);
yunit=90*sin(th);
% 2.3. Inner circles
% inter=[0 1/4.5 1/2.25 1/1.5 1];
inter=[0 1/9 2/9 1/3 4/9 5/9 2/3 7/9 8/9 1];
cax=newplot;
ls=get(cax,'gridlinestyle');
for i=1:length(inter)
   xcircle(i,:)=xunit*inter(i);
   ycircle(i,:)=yunit*inter(i);
end

% 2.4. Filling
fill(xcircle(length(inter),:),ycircle(length(inter),:),'w');
xfill=[90 90 -90 -90 90 90 xcircle(length(inter),:)];
yfill=[0 -90 -90 90 90 0 ycircle(length(inter),:)];
col=get(gcf,'color');
fill(xfill,yfill,col);
fill(xfill,yfill,[1 1 1]);
truc=plot(xfill,yfill,'-');
set(truc,'color',col) ;

% 3. Selection of data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nsel = find(theta>80 | (theta >65 & phi>=270 & phi<=330) ) ;
% V(nsel) = NaN ;

% 4. Projection of the observation directions in pseudo-polar coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = theta.*cos(d2r*phi) ;
y = theta.*sin(d2r*phi) ;

% 5. Min and max of the scale bar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The min and max of the plot are set by the parameters Vmin and Vmax
if Vmax == 0   % then the max of the plot is the max of the data
    nVmax = find(V==max(V(:)));
    Vmax = V(nVmax(1)) ;
end
if Vmin == -1  % then the min of the plot is the min of the data
    nVmin = find(V==min(V(:)));
    Vmin = V(nVmin(1)) ;
end

% 6. Plot with interpolation and extrapolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose between plot_method 1 (recommanded) and plot_method 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_method 1: use "griddata.m" (the surface always goes through the data points)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6.1. Grid (choose between method 1 (default) and 2)

% Grid_method 1 : with frontiers defined by the user
cartsamp = [-90:1:90] ;
[xc,yc] = meshgrid(cartsamp);

% Grid_method 2 : with frontiers from the data
% linx = linspace( min(x(:)) , max(x(:)) ,20) ;
% liny = linspace( min(y(:)) , max(y(:)) ,20) ;
% [xc,yc] = meshgrid(linx,liny);

% 6.2. Resampling

method_griddata = 'natural' ;%changed from cubic  to natural
Vo = griddata(x,y,V,xc,yc,method_griddata) ;

% 6.3. Plot, choose between "contourf" (recommended) and "pcolor"

% Plot with "contourf.m"
nb_contour =10;
contour_values = [Vmin: (Vmax-Vmin)/nb_contour :Vmax] ;
nsel = find(Vo(:) < (Vmin+((Vmax-Vmin)/nb_contour)) ) ; % Just one contour for low values
Vo(nsel) = Vmin ;
contourf(xc,yc,Vo,contour_values);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_method 2: kriging with "geostat_gonio.m" (older version)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6.1. Grid
% cartsamp = [-95:5:95] ;
% [xc,yc] = meshgrid(cartsamp);
%
% 6.2. Resampling "krigeage"
% Vo = geostat_gonio(x,y,V) ;
%
% 6.3. Plot, choose between "contourf" (recommended)and "pcolor"
%
% Plot with "contourf.m"
%%%%%%%%%%%%%%%%%%%%%%%%
% nb_contour = 10;
% contour_values = [Vmin: (Vmax-Vmin)/nb_contour :Vmax];
% nsel = find(Vo(:) < (Vmin+((Vmax-Vmin)/nb_contour)) );
% Vo(nsel) = Vmin;
% contourf(xc,yc,Vo,contour_values)
% [A,B,C] = contourf(xc,yc,Vo,contour_values);
% set(B,'Edgecolor','none');
%
% Plot with "pcolor.m"
%%%%%%%%%%%%%%%%%%%%%%
% h = pcolor(cartsamp,cartsamp,Vo);
% set(h,'EdgeColor','none')
% shading interp

% 7. Add optional features to the plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.1. Colormap
 colormap(1-copper);
% 7.2. Scale bar
caxis([Vmin Vmax]) ;
H_cb = colorbar('vert');
% set(H_cb,'FontSize',16);
% set(H_cb,'FontWeight','bold');
set(H_cb,'FontSize',16,'FontWeight','bold','FontName','Times New Roman','Position',...
    [0.8379,0.122,0.0476,0.6798]);% (changed) Adjust the position of the scale bar to avoid overlapping with the polar plot.
% 7.3. Draw a star for the illumination direction
 plot(thetai*cos(d2r*phii),thetai*sin(d2r*phii),'pk','MarkerSize',16,'MarkerFacecolor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
% 7.4. Plot the sampling points
plot(x,y,'.k','MarkerSize',5) ;
%plot(x,y,'o','MarkerSize',4,'MarkerFacecolor',[0 0 0],'MarkerEdgeColor',[0 0 0]);%

% 7.5. Display the wavelength
text(-90,96,['\lambda = ' num2str(lambda) ' nm ' num2str(thetai) 'бу'],'Fontweight','normal','Fontsize',14,'color','red') ;
% 8. Add features independent of the data points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.1. Outer circle, contour
h=plot(xcircle(length(inter),:),ycircle(length(inter),:),...
    ls,'Color','k','MarkerSize',1);
set(h,'linestyle','-')
% 8.2. Inner circles, contour
for i=2:length(inter)-1
    h=plot(xcircle(i,:),ycircle(i,:),ls,'Color','k','MarkerSize',1);    
end
% 8.3. Phi lines(changed)
for num_axe=1:12
    xmax_axe=cos((num_axe-1)*pi/6)*90; %30 interval
    ymax_axe=sin((num_axe-1)*30*pi/180)*90;
    xaxe_axe(num_axe,:)=linspace(0,xmax_axe,10);
%      xmax_axe=cos((num_axe-1)*pi/4)*90; %45 interval
%   ymax_axe=sin((num_axe-1)*45*pi/180)*90;
%   xaxe_axe(num_axe,:)=linspace(0,xmax_axe,10);
      if (num_axe==1)
        yaxe_axe(num_axe,:)=zeros(1,10);
      else
        yaxe_axe(num_axe,:)=linspace(0,ymax_axe,10);
      end
      plot(xaxe_axe(num_axe,:),yaxe_axe(num_axe,:),ls,'Color','k','MarkerSize',1);
    angle_txt=[int2str((num_axe-1)*30) 'бу'];%
    text(1.1*xmax_axe,1.1*ymax_axe,angle_txt,'horizontalalignment','center','Fontsize',16,'FontName','Times New Roman','FontWeight','bold')
end
% 8.5 (New) Label the viewing zenith angle at 10бу interval
ticks=10;% 
v = [-90 0 0 90]; %
rmin = 0; rmax = v(4); rticks = max(ticks-1,2);
c82 = cos(270*pi/180);% location of label
s82 = sin(270*pi/180);% location of label
rinc = (rmax-rmin)/rticks;
for i=(rmin+rinc):rinc:rmax
    text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
        ['  ' num2str(i)],'verticalalignment','bottom',...
        'handlevisibility','off','FontWeight','bold','Color','k','FontName','Times New Roman','fontsize',12)
end
hold off
set(gcf,'Position',[300 300 550 500]); 
set(gca,'Position',[0.07 0.1 0.72 0.75]);
axis('square')
axis off
hold off