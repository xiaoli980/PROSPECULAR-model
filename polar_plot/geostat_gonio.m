function Vo=geostat_gonio(x,y,V)

% Authors: Stéphane Jacquemoud
% Created: 1998
% Last modified: September 3, 2013
% Affiliation: Institut de Physique du Globe de Paris - University Paris Diderot
%
% Spatial interpolation using the Kriging method

n=length(x);
x=x(:);
y=y(:);
V=V(:);

% calculation of all the distances between the measurements

z=x+y*sqrt(-1);
A=ones(n)*diag(z);
H=abs(A-A.');

% calculation of the covariance function and the variogram

dh=linspace(0,max(H(:))/2,5);
delta=(dh(2)-dh(1))/2;
k=[];
c=ones(1,5);
g=ones(1,5);
for j=1:5
  if dh(j)==0
    [Nx,Ny]=find(H<=delta);
  else
    [Nx,Ny]=find((H>dh(j)-delta & H<=dh(j)+delta));
  end
  if length(Nx)==0
    k=[k j];
  else
    c(j)=mean((V(Nx)-mean(V(Nx))).*(V(Ny)-mean(V(Ny))));
    g(j)=mean((V(Nx)-V(Ny)).^2)/2;
  end
end
dh(k)=[];
c(k)=[];
g(k)=[];

% fitting of the variogram with a straight line

lin=polyfit(dh,g,1);
  
% estimation of the variables on a 39-by-39 regular grid

C=[lin(1)*H+lin(2) -ones(n,1) ; -ones(1,n)  0];
c_x=1;

for xo=-95:5:95
  c_y=1;
  for yo=-95:5:95
    zo=xo+yo*sqrt(-1);
    Ho=abs(z-zo);
    D=[lin(1)+lin(2)*Ho ; -1];
    wo=C\D;
    Vo(c_y,c_x)=sum(wo(1:n).*V);
    c_y=c_y+1;
  end
  c_x=c_x+1;
end