azimuth = 352;     %here to set Sigma SOUTH
sigma1 = -50;
sigma2 = -20;
sigma3 = -20;
az = azimuth-45; %transfer to simulation coordinate by rotation of 45 degree
sigma_yy_SOUTH = sigma1 * cosd(az) *cosd(az)+ sigma2 * sind(az) * sind(az);
sigma_xx_SOUTH = sigma2 * cosd(az) *cosd(az)+ sigma1 * sind(az) * sind(az);
sigma_xy_SOUTH = (sigma1-sigma2) * cosd(az) * sind(az);
sigma_zz_SOUTH = sigma3;
azimuth = 320;    %here to set Sigma North
sigma1 = -50;
sigma2 = -20;
sigma3 = -20;
az=azimuth-45;
sigma_yy_NORTH = sigma1 * cosd(az) *cosd(az)+ sigma2 * sind(az) * sind(az);
sigma_xx_NORTH = sigma2 * cosd(az) *cosd(az)+ sigma1 * sind(az) * sind(az);
sigma_xy_NORTH = (sigma1-sigma2) * cosd(az) * sind(az);
sigma_zz_NORTH = sigma3;

sigma_xx = sigma_xx_NORTH * min((yy+150)/(200),1) + sigma_xx_SOUTH * max((50-yy)/200,0);
sigma_yy = sigma_yy_NORTH * min((yy+150)/(200),1) + sigma_yy_SOUTH * max((50-yy)/200,0);
sigma_xy = sigma_xy_NORTH * min((yy+150)/(200),1) + sigma_xy_SOUTH * max((50-yy)/200,0);
sigma_zz = sigma_zz_NORTH * min((yy+150)/(200),1) + sigma_zz_SOUTH * max((50-yy)/200,0);


data=FSEM3D_snapshot(0);
[yy,zz]=meshgrid(min(data.Y):0.1:max(data.Y),min(data.Z):0.1:max(data.Z));
xx=griddata(data.Y,data.Z,data.X,yy,zz);

[dxdy,dxdz] = gradient(xx,0.1);
norm=sqrt(dxdy.^2+dxdz.^2+1);
ny = -dxdy./norm;
nz = -dxdz./norm;
nx = 1./norm;
sx = -ny;
sy = nx;
sz = zeros(size(sx));

norms = sqrt(sx.^2+sy.^2);
sx = sx ./ norms;
sy = sy ./ norms;

dx = -sy .* nz;
dy = sx .* nz;
dz = sy .* nx - sx .* ny;
normd = sqrt(dx.^2+dy.^2 +dz.^2);
dx = dx ./normd;
dy = dy ./normd;
dz = dz ./normd;

ss = sigma_xx.*nx.*sx + sigma_xy.*ny.*sx+...
     sigma_xy.*nx.*sy + sigma_yy.*ny.*sy+sigma_zz.*nz.*sz;    %compute the along strike shear stress
ds = sigma_xx.*nx.*dx + sigma_xy.*ny.*dx+...
     sigma_xy.*nx.*dy + sigma_yy.*ny.*dy+sigma_zz.*nz.*dz;    %compute the slong dip shear stress
ns = sigma_xx.*nx.*nx + sigma_xy.*ny.*nx+...
     sigma_xy.*nx.*ny + sigma_yy.*ny.*ny+sigma_zz.*nz.*nz;     %compute normal stress
 


