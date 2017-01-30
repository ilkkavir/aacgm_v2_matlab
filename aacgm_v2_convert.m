function [lat,lon,r] = aacgm_v2_convert(in_lat,in_lon,in_height,time,code,coord)
%
% [lat,lon,r] = aacgm_v2_convert(in_lat,in_lon,in_height,time,code,coord)
% 
% Conversion between geographic and altitude-adjusted
% geomagnetic coordinates. This function and everything this calls
% are modified from the original C version by Shepherd. See
% Shepherd, S.G., "Altitude-Adjusted Corrected Geomagnetic
% Coordinates: Definition and Functional Approximations", JGR, 119,
% 7501-7521, 2014. DOI: 10.1002/2014JA020264.
% 
% See the function convert_geo_coord for details
%
% INPUT:
%   in_lat     latitude (deg)
%   in_lon     longitude (deg)
%   in_height  height(km)
%   time       time as MATLAB datetime object
%   code       0 for (geographic) to AACGM-v2, 1 for the inverse
%              transform
%   coord      geographic coordinate system, 0 for geodetic, 1 for
%              geocentric. THE MAGNETIC COORDINATES ARE ALWAYS
%              GEOCENTRIC! coord defines the input coordinate
%              system when code=0 and the output system when code=1.
%
% 
% OUTPUT:
%   lat    latitude (deg)
%   lon    longitude (deg)
%   r      height (km)
% 
% IV 2016
% 

persistent RE MAXALT
if isempty(RE) | isempty(MAXALT)
    load('aacgmv2coefs.mat','RE','MAXALT')
end

if in_height < 0
    error(['coordinate transformations are not intended for altitudes ' ...
           '< 0 km: ',num2str(in_height)])
end
if in_height > MAXALT
    error(['coordinate transformations are not intended for altitudes ' ...
           '> ',num2str(MAXALT),' km: ',num2str(in_height)])
end

% /* latitude out of bounds */
if abs(in_lat) > 90
    error(['latitude must be in the range of -90 to 90 degrees: ',num2str(in_lat)]);
end


% if forward calculation, we must convert the input coordinates
% into geocentric if they are geodetic
if code==0 & coord==0
    [lati,long,heig] = geod2geoc(in_lat,in_lon,in_height,RE);
else
    lati = in_lat;
    long = in_lon;
    heig = in_height;
end

% time as decimal year
dyear = (posixtime(time)-posixtime(datetime(year(time),1,1)))/ ...
        (posixtime(datetime(year(time)+1,1,1))- ...
         posixtime(datetime(year(time),1,1))) + year(time);

% the actual conversion
[lat_out,lon_out] = convert_geo_coord(lati,long,heig,dyear,code);

% coordinate conversion if necessary
if code==1 & coord==0
    [lat_out,lon_out,r] = geoc2geod(lat_out,lon_out,heig,RE);
else 
    r = heig;
end

lat = lat_out;
lon = lon_out;
r = r;

end 









function [lat_out,lon_out] = convert_geo_coord(lat_in,lon_in, ...
                                               height_in,dyear,code)
% 
% [lat_out,lon_out] =
% convert_geo_coord(lat_in,lon_in,height_in,dyear,code) 
%
% Conversion between geocentric and altitude-adjusted
% geomagnetic coordinates. This function and everything this calls
% are modified from the original C version by Shepherd. See
% Shepherd, S.G., "Altitude-Adjusted Corrected Geomagnetic
% Coordinates: Definition and Functional Approximations", JGR, 119,
% 7501-7521, 2014. DOI: 10.1002/2014JA020264.
% 
% The tabulated spherical harmonic coefficients and the original
% source code were downloaded from
% https://engineering.dartmouth.edu/superdarn/aacgm.html 
%
% This is a significantly reduced version of the conversions
% introduced in the paper. Only the spherical harmonics valid below
% 2000 km in altitude are implemented.
%
% This function reads the spherical harmonic coefficients from file 
% aacgmv2coefs.mat. The file can be updated with the function
% update_aacgmv2coefs.  
%
% INPUT:
%   lat_in     latitude (deg)
%   lon_in     longitude (deg)
%   height_in  height(km)
%   dyear      time as decimal year (e.g. 2014.5)
%   code       0 for (geodetic) to AACGM-v2, 1 for the inverse transform
% 
% OUTPUT:
%   lat_out    latitude (deg)
%   lon_out    longitude (deg)
% 
% IV 2016
% 

DTOR = pi/180; % degrees to radians...

persistent aacgmv2coefs aacgmv2years MAXALT AACGM_KMAX RE SHORDER

% load the tabulated coefficients if they were not yet loaded
if isempty(aacgmv2coefs) | isempty(aacgmv2years) | isempty(MAXALT) | isempty(AACGM_KMAX) | isempty(RE) | isempty(SHORDER)
    load('aacgmv2coefs.mat');
end

if dyear > max(aacgmv2years) | dyear < min(aacgmv2years)
    error(['dyear must be between ',num2str(min(aacgmv2years)),' and ', num2str(max(aacgmv2years))]);
end

% normalized altitude and its powers
alt_var = height_in / MAXALT;
alt_var_sq = alt_var^2;
alt_var_cu = alt_var * alt_var_sq;
alt_var_qu = alt_var * alt_var_cu;

% interpolate the spherical harmonic coefficients
% linear interpolation in time

[yytmp,iitmp] = sort(abs(aacgmv2years-dyear));
iiprev = min(iitmp(1:2));
iinext = max(iitmp(1:2));
yyprev = aacgmv2years(iiprev);
yynext = aacgmv2years(iinext);
weightprev = (dyear-yyprev)/(yynext-yyprev);
weightnext = 1 - weightprev;
coefsdyear = squeeze(aacgmv2coefs(iiprev,:,:,:,:).*weightprev + ...
                     aacgmv2coefs(iinext,:,:,:,:).*weightnext);

% quadric polynomial in height
cint = coefsdyear(:,:,1,:) + coefsdyear(:,:,2,:).*alt_var + ...
       coefsdyear(:,:,3,:).*alt_var_sq + coefsdyear(:,:,4,:).*alt_var_cu ...
       + coefsdyear(:,:,5,:).*alt_var_qu;


% the actual transform begins here
x = 0;
y = 0;
z = 0;

lon_input = lon_in*DTOR;

if code==1 % inverse transform aacgmv2 -> geodetic
    r1 = cos(lat_in*DTOR);
    ra = (height_in/RE + 1)*(r1*r1);
    if (ra > 1)
        error()
    end
    r1 = acos(sqrt(ra));
    lat_adj = abs(r1)*sign(lat_in);
end
colat_input = (90.-lat_in)*DTOR;


% compute the values of spherical harmonic functions. There seems
% to be some non-standard normalization involved. To be sure that
% this will be OK,  I have translated the original C-function to
% MATLAB. IV.
ylmval = AACGM_v2_Rylm(colat_input, lon_input, SHORDER , AACGM_KMAX);
for l=0:SHORDER
    for m=-l:l
        k = l * (l+1) + m + 1;
        x = x + cint(k,1,code+1)*ylmval(k);
        y = y + cint(k,2,code+1)*ylmval(k);
        z = z + cint(k,3,code+1)*ylmval(k);
    end
end



if code == 0
    fac = x*x + y*y;
    if fac > 1.
        % /* we are in the forbidden region and the solution is undefined */
        lat_out = NaN;
        lon_out = NaN;
        return
    end
    
    ztmp = sqrt(1. - fac);
    z = sign(z)*ztmp;

else 

    %/* SGS - for inverse the old normalization produces lower overall errors...*/
    r = sqrt(x*x + y*y + z*z);
    if ((r< 0.9) | (r > 1.1))
        lat_out = NaN;
        lon_out = NaN;
        return
    end
    z = z/r;
    x = x/r;
    y = y/r;
end 

if (z > 1.)
    colat_temp = 0;
elseif (z < -1.)
    colat_temp = M_PI;
else
    colat_temp = acos(z);
end

if ((abs(x) < 1e-8) && (abs(y) < 1e-8))
    lon_temp = 0;
else 
    lon_temp = atan2(y,x);
end

lon_output = lon_temp;

colat_output = colat_temp;

lat_out = 90. - colat_output/DTOR;
lon_out = lon_output/DTOR;

end 



















function ylmval = AACGM_v2_Rylm( colat, lon, order, nylmval)
%
% MATLAB version of the C-function AACGM_v2_Rylm.
%
% I do not exactly understand everything here, so I tried to avoid
% changing anything, the result does not look nice but works...
% 
% IV 2016
%

ylmval = zeros(nylmval,1);

cos_theta = cos(colat);
sin_theta = sin(colat);

cos_lon = cos(lon);
sin_lon = sin(lon);

d1 = -sin_theta;
z2x = cos_lon;
z2y = sin_lon;

z1x = d1 * z2x;
z1y = d1 * z2y;
q_facx = z1x;
q_facy = z1y;

%    /*
%     * Generate Zonal Harmonics (P_l^(m=0) for l = 1,order) using recursion
%     * relation (6.8.7), p. 252, Numerical Recipes in C, 2nd. ed., Press. W.
%     * et al. Cambridge University Press, 1992) for case where m = 0.
%     *
%     * l Pl = cos(theta) (2l-1) Pl-1 - (l-1) Pl-2 (6.8.7)
%     *
%     * where Pl = P_l^(m=0) are the associated Legendre polynomials
%     *
%     */

ylmval(1) = 1; %/* l = 0, m = 0 */
ylmval(3) = cos_theta; %/* l = 1, m = 0 */

for l=2:order
    %/* indices for previous two values: k = l * (l+1) + m with m=0 */
    ia = (l-2)*(l-1) + 1;
    ib = (l-1)*l + 1;
    ic = l * (l+1) + 1;
    
    ylmval(ic) = (cos_theta * (2*l-1) * ylmval(ib) - (l-1)* ...
                  ylmval(ia))/l;
end
% /*
%  * Generate P_l^l for l = 1 to (order+1)^2 using algorithm based upon (6.8.8)
%  * in Press et al., but incorporate longitude dependence, i.e., sin/cos (phi)
%  *                              
%  *
%  * Pll = (-1)^l (2l-1)!! (sin^2(theta))^(l/2)
%  *
%  * where Plm = P_l^m are the associated Legendre polynomials
%  *
% */
q_valx = q_facx;
q_valy = q_facy;
ylmval(4) = q_valx; %/* l = 1, m = +1 */
ylmval(2) = -q_valy; %/* l = 1, m = -1 */
for l=2:order
    d1 = l*2 - 1.;
    z2x  = d1 * q_facx;
    z2y  = d1 * q_facy;
    z1x  = z2x * q_valx - z2y * q_valy;
    z1y  = z2x * q_valy + z2y * q_valx;
    q_valx = z1x;
    q_valy = z1y;
    
    % /* indices for previous two values: k = l * (l+1) + m */
    ia = l*(l+2) + 1; %/* m = +l */
    ib = l*l + 1;     %/* m = -l */
    ylmval(ia) =  q_valx;
    ylmval(ib) = -q_valy;
end

%/*
% * Generate P_l,l-1 to P_(order+1)^2,l-1 using algorithm based upon (6.8.9)
% * in Press et al., but incorporate longitude dependence, i.e., sin/cos (phi)
% *
% * Pl,l-1 = cos(theta) (2l-1) Pl-1, l-1
% *
% */
 
for l=2:order
    l2 = l*l;
    tl = 2*l;
    %    /* indices for Pl, l-1; Pl-1,l-1; Pl,-(l-1); Pl-1,-(l-1) */
    ia = l2 - 1 + 1;
    ib = l2 - tl + 1 + 1;
    ic = l2 + tl - 1 + 1;
    id = l2 + 1 + 1;
    
    fac = tl - 1;
    ylmval(ic) = fac * cos_theta * ylmval(ia); %/* Pl,l-1   */
    ylmval(id) = fac * cos_theta * ylmval(ib); %/* Pl,-(l-1) */
end

% /*
% * Generate remaining P_l+2,m to P_(order+1)^2,m for each m = 1 to order-2
% * using algorithm based upon (6.8.7) in Press et al., but incorporate
% * longitude dependence, i.e., sin/cos (phi).
% *
% * for each m value 1 to order-2 we have P_mm and P_m+1,m so we can compute
% * P_m+2,m; P_m+3,m; etc.
% *
% */

for m=1:order-2
    for l=m+2:order
        ca = (2*l-1)/(l-m);
        cb = (l+m-1)/(l-m);
        
        l2 = l*l;
        ic = l2 + l + m + 1;
        ib = l2 - l + m + 1;
        ia = l2 - l - l - l + 2 + m + 1;
        % /* positive m */
        ylmval(ic) = ca * cos_theta * ylmval(ib) - cb * ylmval(ia);
        
        ic = ic - (m+m);
        ib = ib - (m+m);
        ia = ia - (m+m);
        % /* negative m */
        ylmval(ic) = ca * cos_theta * ylmval(ib) - cb * ylmval(ia);
    end
end

%/*
% * Normalization added here (SGS)
% * 
% * Note that this is NOT the standard spherical harmonic normalization factors
% * 
% * The recursive algorithms above treat positive and negative values of m in
% * the same manner. In order to use these algorithms the normalization must
% * also be modified to reflect the symmetry.
% * 
% * Output values have been checked against those obtained using the internal
% * IDL legendre() function to obtain the various associated legendre
% * polynomials.
% * 
% * As stated above, I think that this normalization may be unnecessary. The
% * important thing is that the various spherical harmonics are orthogonal,
% * rather than orthonormal.
% * 
% */

% /* determine array of factorials */
% /*fact = (long unsigned *)malloc(sizeof(long unsigned)*(2*order+2));*/
fact = zeros((2*order+2),1);
fact(1) = 1;
fact(2) = 1;
for k=2:(2*order+1)
    fact(k+1) = k*fact(k); 
end

ffff = zeros((order+1)*(order+1),1);

% /* determine normalization factors */
for l=0:order
    for m=0:l
        k = l * (l+1) + m; %/* 1D index for l,m */
        ffff(k+1) = sqrt((2*l+1)/(4*pi) * fact(l-m+1)/fact(l+m+1));
        ylmval(k+1) = ylmval(k+1) * ffff(k+1);
    end
    for m=-l:-1
        k  = l * (l+1) + m; %/* 1D index for l,m */
        kk = l * (l+1) - m;
        ylmval(k+1) = ylmval(k+1) * ffff(kk+1); %* ((-m % 2) ? -1
                                                %      % : 1);
        if mod(-m,2)~=0
            ylmval(k+1) = -1*ylmval(k+1);
        end
    end
end

end




function [lat_out,lon_out,r_out] = geoc2geod(lat,lon,r,RE)
%
% [lat_out,lon_out,r_out] = geoc2geoe(lat,lon,r,RE)
%
%  Conversion from geocentric to geodetic coordinates. Modified
%  from the c-function geoc2geoe in igrflib.c
%
% INPUT:
%  lat geocentric latitude (deg)
%  lon geocentric longitude (deg)
%  alt geocentric altitude (radial distance - RE) (km)
%  RE  Earh radius (km)
%
% OUTPUT:
%  lat_out geodetic latitude (deg)
%  lon_out geodetic longitude (deg)
%  r_out   geodetic altitude (km)
%
% IV 2016
%

DTOR = pi/180;
a = 6378.1370; %            /* semi-major axis */
f = 1./298.257223563; %    /* flattening */
b = a*(1. -f);       %      /* semi-minor axis */
ee = (2. - f) * f;
e4 = ee*ee;
aa = a*a;

theta = (90. - lat)*DTOR;
phi   = lon * DTOR;

st = sin(theta);
ct = cos(theta);
sp = sin(phi);
cp = cos(phi);

x = (r+RE) * st * cp;
y = (r+RE) * st * sp;
z = (r+RE) * ct;

k0i   = 1. - ee;
pp    = x*x + y*y;
zeta  = k0i*z*z/aa;
rho   = (pp/aa + zeta - e4)/6.;
s     = e4*zeta*pp/(4.*aa);
rho3  = rho*rho*rho;
t     = power(rho3 + s + sqrt(s*(s+2*rho3)), 1./3.);
u     = rho + t + rho*rho/t;
v     = sqrt(u*u + e4*zeta);
w     = ee*(u + v - zeta)/(2.*v);
kappa = 1. + ee*(sqrt(u+v+w*w) + w)/(u + v);

lat_out = atan2(z*kappa,sqrt(pp))/DTOR;
lon_out = lon;
r_out = sqrt(pp + z*z*kappa*kappa)/ee * (1./kappa - k0i);


end % geoc2geod





function [lat_out,lon_out,r_out] = geod2geoc(lat,lon,alt,RE)
%
% [lat_out,lon_out,r_out] = geod2geoc(lat,lon,alt,RE)
%
%  Conversion from geodetic to geocentric coordinates. Modified
%  from the c-function geod2geoc in igrflib.c
%
% INPUT:
%  lat geodetic latitude (deg)
%  lon geodetic longitude (deg)
%  alt geodetic altitude (km)
%  RE  Earh radius (km)
%
% OUTPUT:
%  lat_out geocentric latitude (deg)
%  lon_out geocentric longitude (deg)
%  r_out   geocentric altitude (radial distance - RE) (km)
%
% IV 2016
%

DTOR = pi/180;

a = 6378.1370;%/* semi-major axis */
f = 1./298.257223563; %/* flattening */
b = a*(1. -f);%/* semi-minor axis */
a2 = a*a;
b2 = b*b;
theta = (90. -lat)*DTOR;%/* colatitude in radians   */
st = sin(theta);
ct = cos(theta);
one = a2*st*st;
two = b2*ct*ct;
three = one + two;
rho = sqrt(three);%/* [km] */
r = sqrt(alt*(alt+2*rho) + (a2*one + b2*two)/three);%    /* [km] ...*/
cd = (alt+rho)/r;
sd = (a2-b2)/rho *ct*st/r;

r_out = r - RE;
lat_out = 90 - acos(ct*cd - st*sd)/DTOR;
lon_out = lon;

end % geod2geoc
