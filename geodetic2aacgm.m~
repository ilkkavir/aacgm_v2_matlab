function [lat,lon,r] = geodetic2aacgm(inlat,inlon,inheight,time)
% 
% 
% [lat,lon,r] = geodetic2aacgm(inlat,inlon,inheight,time)
% 
% Conversion from geodetic (WGS84) to altitude adjusted geomagnetic
% coordinates. This is a wrapper to aacgm_v2_convert.
% 
% INPUT:
%  inlat     geocentric latitude (deg)
%  inlon     geocentric longitude (deg)
%  inheight  ellipsoid height [km]
%  time      time as matlab datetime structure
% 
% OUTPUT: 
%  lat       geomagnetic latitude (deg)
%  lon       geomagnetic longitude (deg)
%  r         geocentric height (radial distance - RE ) (km)
%
% 
% See alo geocentric2aacgm, aacgm2geocentric, aacgm2geodetic, aacgm_v2_convert
% 
% IV 2016
%

[lat lon r] = aacgm_v2_convert(inlat,inlon,inheight,time,0,0);

end