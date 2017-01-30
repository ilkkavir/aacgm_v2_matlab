function [mlt,slat,slon] = ut2mlt_mlon(time,mlon)
% 
% 
%  [mlt,slat,slon] = ut2mlt_mlon(time,mlon)
% 
% Conversion from UT time to magnetic local time. This is a wrapper
% to magneticLocalTime.
%
% INPUT:
%  time   time as MATLAB datetime
%  mlon    magnetic longitude (deg)
%
% OUTPUT:
%  mlt    magnetic local time (h)
%  slat   latitude of the subsolar point in geocentric coordinates (deg)
%  slon   longitude of the subsolar point in geocentric coordinates (deg)
%
% 
% Seel alo ut2mlt_geodetic, ut2mlt_geocentric
% 
% IV 2016
%


[mlt,slat,slon] = magneticLocalTime(time,mlon);

end