function [lat,lon,r] = aacgm2geodetic(inlat,inlon,inheight,time)
% 
% 
% [lat,lon,r] = aacgm2geodetic(inlat,inlon,inheight,time)
% 
% Conversion from altitude adjusted geomagnetic to geodetic (WGS84)
% coordinates. This is a wrapper to aacgm_v2_convert.
% 
% INPUT:
%  inlat     geomagnetic latitude (deg)
%  inlon     geomagnetic longitude (deg)
%  inheight  geocentric height (radial distance - RE ) (km)
%  time      time as matlab datetime structure
% 
% OUTPUT: 
%  lat       geodetic latitude (deg)
%  lon       geodetic longitude (deg)
%  r         ellipsoid height [km]
%
% 
% See alo geodetic2aacgm, geocentric2aacgm, aacgm2geocentric, aacgm_v2_convert
% 
% IV 2016
%

[lat lon r] = aacgm_v2_convert(inlat,inlon,inheight,time,1,0);

end