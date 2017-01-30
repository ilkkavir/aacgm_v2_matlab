function [lat,lon,r] = aacgm2geocentric(inlat,inlon,inheight,time)
% 
% 
% [lat,lon,r] = aacgm2geocentric(inlat,inlon,inheight,time)
% 
% Conversion from altitude adjusted geomagnetic to geocentric
% coordinates. This is a wrapper to aacgm_v2_convert.
% 
% INPUT:
%  inlat     geomagnetic latitude (deg)
%  inlon     geomagnetic longitude (deg)
%  inheight  geocentric height (radial distance - RE ) (km)
%  time      time as matlab datetime structure
% 
% OUTPUT: 
%  lat       geocentric latitude (deg)
%  lon       geocentric longitude (deg)
%  r         geocentric height (radial distance - RE ) (km)
%
% 
% See alo geodetic2aacgm, geocentric2aacgm, aacgm2geodetic, aacgm_v2_convert
% 
% IV 2016
%

[lat lon r] = aacgm_v2_convert(inlat,inlon,inheight,time,1,1);

end