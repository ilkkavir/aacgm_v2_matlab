function [mlt,slat,slon] = ut2mlt_geodetic(time,inlat,inlon,inheight)
% 
% 
%  [mlt,slat,slon] = ut2mlt_geodetic(time,inlat,inlon,inheight)
% 
% Conversion from UT time to magnetic local time, with location
% given in geodetic (WGS84) coordinates.
% 
% INPUT:
%  time     time as matlab datetime struct
%  inlat    geodetic latitude (deg)
%  inlon    geodetic longitude (deg)
%  inheight ellipsoid height (km)
% 
% OUTPUT: 
%  mlt     magnetic local time (hours)
%  slat   latitude of the subsolar point in geocentric coordinates (deg)
%  slon   longitude of the subsolar point in geocentric coordinates (deg)
%
% 
% Seel alo ut2mlt_mlon, ut2mlt_geocentric
% 
% IV 2016
%

[mlat mlon mr] = aacgm_v2_convert(inlat,inlon,inheight,time,0,0);

[mlt,slat,slon] = magneticLocalTime(time,mlon);

end