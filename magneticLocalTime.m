function [mlt,slat,slon] = magneticLocalTime(time,mlon)
%
% mlt = magneticLocalTime(time,lon)
%
% Calculate magnetic local time and geocentric coordinates of the
% subsolar point. Based on MLTConvert_v2 in aacgm_v2.
% Functions called by magneticLocalTime are based on routines in
% AstAlglib, see LICENSE-AstAlg.txt for details. Converted to
% MATLAB by IV.
%
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
% IV 2016
%

DTOR = pi/180;

persistent time_prev mlon_prev mlt_prev slat_prev slon_prev

% do not recalculate if the function was called with identical
% arguments as previously
if ~isempty(time_prev)
    if time_prev == time
        if ~isempty(mlon_prev)
            if mlon_prev == mlon
                if ~isempty(mlt_prev) & ~isempty(slat_prev) & ~isempty(slon_prev)
                    mlt = mlt_prev;
                    slat = slat_prev;
                    slon = slon_prev;
                    return
                end
            end
        end
    end
end


hgt = 700;

% julian day (matlab builtin!)
jd = juliandate(time);

% equation of time
eqt = AstAlg_equation_of_time(jd);

% solar declination (latitude of subsolar point)
dec = AstAlg_solar_declination(jd);

% UT hour
ut = time.Hour*3600 + time.Minute*60 + time.Second;

% correction from eqution of time
at = ut + eqt*60;

% subsolar point longitude
slon = (43200 - at)*15/3600;

% aacgm-v2 coordinates of the reference point
[mlat,mlon_ref,r] = aacgm_v2_convert(dec,slon,hgt,time,0,0);

% mlt based on subsolar point
mlt = 12 + ( mlon - mlon_ref )/15;
mlt = mod(mlt,24);

slat = dec;
time_prev = time;
mlon_prev = mlon;
mlt_prev = mlt;
slat_prev = dec;
slon_prev = slon;

end


function eqt = AstAlg_equation_of_time(jd)
%
% eqt = AstAlg_equation_of_time(jd)
%
% Equation of time, i.e. difference between apparent and mean solar
% time. This may be up to 20 minutes!
%
% Based on the c-function AstAlg_equation_of_time in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  eqt equation of time. A positive value menas that the true sun
%  crosses the observer's meridian before the mean sun (i.e. mean
%  time is lagging).
%
% IV 2016

persistent last_jd DTOR last_eqt

if isempty(DTOR)
    DTOR=pi/180;
end

% if the value was already calculated, simply return it
if ~isempty(last_jd)
    if jd==last_jd;
        if ~isempty(last_eqt)
            eqt = last_eqt;
            return
        end
    end
end

% mean solar longitude
sml = AstAlg_mean_solar_longitude(jd);
% solar right ascension
sra = AstAlg_solar_right_ascension(jd);
% mean obliquity of the Earth
obliq = AstAlg_mean_obliquity(jd);
% nutation correction
[dpsi,deps] = AstAlg_nutation_corr(jd);


% put these together
eqt = sml - 0.0057183 - sra + dpsi * cos( DTOR*(obliq + deps) );

eqt = mod(eqt,360);

% from degrees to minutes
eqt = 4 * eqt;

if eqt > 20.0
    eqt = eqt - 24*60; % wrap back 24 hours
end
if eqt < -20.0
    eqt = eqt + 24*60; % wrap forward 24 hours
end

last_jd = jd;
last_eqt = eqt;

end



function sml = AstAlg_mean_solar_longitude(jd)
%
% sml = AstAlg_mean_solar_longitude(jd)
%
% Mean solar longitude.
%
% Based on the c-function AstAlg_mean_solar_longitude in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  sml mean solar longitude (deg)
%
% IV 2016

persistent last_jd coefs last_sml

if isempty(coefs)
    coefs = [280.4664567, 360007.6982779,0.03032028, 2.00276381406e-5, ...
             -6.53594771242e-5, -0.50e-6];
end

if ~isempty(last_jd)
    if last_jd == jd
        if ~isempty(last_sml)
            sml = last_sml;
            return
        end
    end
end


tau = ( jd - J2000() )/365250;

% solar longitude
sml = 0;
for i=6:-1:1
    sml = tau * sml + coefs(i);
end

sml = mod(sml,360);

last_jd = jd;
last_sml = sml;

end

function jd = J2000()
%
% Julian date for the J2000 epoch
%
jd = 2451545;

end



function sra = AstAlg_solar_right_ascension(jd);
%
% sra = AstAlg_solar_right_ascension(jd);
%
% Solar right ascension
%
% Based on the c-function AstAlg_solar_right_ascension in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  sra solar right ascension (deg)
%
% IV 2016

persistent last_jd last_sra DTOR
if isempty(DTOR)
    DTOR = pi/180;
end

if ~isempty(last_jd)
    if last_jd == jd
        if ~isempty(last_sra)
            sra = last_sra;
            return
        end
    end
end

% solar longitude in radians
slong = DTOR * AstAlg_apparent_solar_longitude(jd);

% obliquity in radians
epsil = DTOR * AstAlg_apparent_obliquity(jd);

alpha = atan2(cos(epsil)*sin(slong),cos(slong));

sra = alpha/DTOR;

last_jd = jd;
last_sra = sra;

end


function asl = AstAlg_apparent_solar_longitude(jd);
%
% asl = AstAlg_apparent_solar_longitude(jd);
%
% Apparent solar longitude
%
% Based on the c-function AstAlg_apparent_solar_longitude in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  asl apparent solar longitude (deg)
%
% IV 2016
persistent last_jd last_asl DTOR
if isempty(DTOR)
    DTOR = pi/180;
end

if ~isempty(last_jd)
    if last_jd == jd
        if ~isempty(last_asl)
            asl = last_asl;
            return
        end
    end
end

gsl = AstAlg_geometric_solar_longitude(jd);
lan = AstAlg_lunar_ascending_node(jd);

asl = gsl - 0.00569 - 0.00478 * sin( DTOR * lan );

last_jd = jd;
last_asl = asl;

end




function gsl = AstAlg_geometric_solar_longitude(jd);
%
% gsl = AstAlg_geometric_solar_longitude(jd);
%
% Geometric solar longitude
%
% Based on the c-function AstAlg_geometric_solar_longitude in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  gsl geometric solar longitude (deg)
%
% IV 2016
persistent last_jd last_gsl DTOR
if isempty(DTOR)
    DTOR = pi/180;
end

if ~isempty(last_jd)
    if last_jd == jd
        if ~isempty(last_gsl)
            gsl = last_gsl;
            return
        end
    end
end

% time difference since J2000 in centuries...
tau = ( jd - J2000() ) / 36525;

% solar mean longitude
sml = AstAlg_mean_solar_longitude(jd);

% mean solar anomaly
sma = AstAlg_mean_solar_anomaly(jd);

smar = DTOR * sma;
gc = (1.914602 - 0.004817*tau - 0.000014*(tau*tau)) * sin(smar) + ...
     (0.019993 - 0.000101*tau) * sin(2.0*smar) + 0.000289 * sin(3.0*smar);

gsl = sml + gc;

gsl = mod(gsl,360);

last_jd = jd;
last_gsl = gsl;

end





function sma = AstAlg_mean_solar_anomaly(jd);
%
% sma = AstAlg_mean_solar_anomaly(jd);
%
% Mean solar anomaly
%
% Based on the c-function AstAlg_mean_solar_anomaly in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  sma mean solar anomaly (deg)
%
% IV 2016
persistent last_jd last_sma

if ~isempty(last_jd)
    if last_jd == jd
        if ~isempty(last_sma)
            sma = last_sma;
            return
        end
    end
end

tau = ( jd - J2000() ) / 36525;

sma = 357.5291130 + 35999.05029 * tau - 0.0001537 * (tau*tau);

sma = mod(sma,360);

last_jd = jd;
last_sma = sma;

end




function lan = AstAlg_lunar_ascending_node(jd);
%
% lan = AstAlg_lunar_ascending_node(jd);
%
% Lunar ascending node
%
% Based on the c-function AstAlg_lunar_ascending_node in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  lan lunar ascending node (deg)
%
% IV 2016
persistent last_jd last_lan

if ~isempty(last_jd)
    if last_jd == jd
        if ~isempty(last_lan)
            lan = last_lan;
            return
        end
    end
end

% time difference from J2000 in centuries..
tau = ( jd - J2000() )  / 36525;

omega = (((tau/4.50e5 + 2.0708e-3)*tau - 1.934136261e3)*tau) + ...
        125.04452;

lan = mod(omega,360);

last_jd = jd;
last_lan = lan;

end



function aob = AstAlg_apparent_obliquity(jd);
%
% aob = AstAlg_apparent_obliquity(jd);
%
% Apparent obliquity (angle between Earth's spin axis and the ecliptic)
%
% Based on the c-function AstAlg_apparent_obliquity in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  aob apparent obliquity (deg)
%
% IV 2016
persistent last_jd last_aob DTOR
if isempty(DTOR)
    DTOR = pi/180;
end

if ~isempty(last_jd)
    if last_jd == jd
        if ~isempty(last_aob)
            aob = last_aob;
            return
        end
    end
end

mob = AstAlg_mean_obliquity(jd);
lan = AstAlg_lunar_ascending_node(jd);

aob = mob + 0.00256 * cos( DTOR * lan );

last_jd = jd;
last_aob = aob;

end




function mob = AstAlg_mean_obliquity(jd);
%
% mob = AstAlg_mean_obliquity(jd);
%
% Mean obliquity (angle between Earth's spin axis and the ecliptic)
%
% Based on the c-function AstAlg_mean_obliquity in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  mob mean obliquity (deg)
%
% IV 2016
persistent last_jd last_mob coefs
if isempty(coefs)
    coefs = [23.439291111111, -0.0130041666667, -1.638888889e-7, 5.036111111e-7];
end

if ~isempty(last_jd)
    if last_jd == jd
        if ~isempty(last_mob)
            mob = last_mob;
            return
        end
    end
end

tau = ( jd - J2000() ) / 36525;

mob = ((((coefs(4)*tau) + coefs(3)) * tau) + coefs(2)) * tau + coefs(1);

last_jd = jd;
last_mob = mob;

end



function [slong_corr,obliq_corr] = AstAlg_nutation_corr(jd)
%
% [slong_corr,obliq_corr] = AstAlg_nutation_corr(jd)
%
% Correction to solar longitude and obliquity due to nutation
%
% Based on the c-function AstAlg_nutation_corr in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  slong_corr correction to solar longitude (deg)
%  obliq_corr correction to obliquity (deg)
%
% IV 2016

persistent last_jd last_slong_corr last_obliq_corr DTOR
if isempty(DTOR)
    DTOR = pi/180;
end
if ~isempty(last_jd)
    if last_jd == jd
        if ~isempty(last_slon_corr) & ~isempty(last_obliq_corr)
            slon_corr = last_slon_corr;
            obliq_corr = last_obliq_corr;
            return
        end
    end
end

% mean solar longitude
msl = AstAlg_mean_solar_longitude(jd);

% lunar longitude
mll = AstAlg_mean_lunar_longitude(jd);

% lunar ascending node
lan = AstAlg_lunar_ascending_node(jd);

% in radians
mslr = DTOR*msl;
mllr = DTOR*mll;
lanr = DTOR*lan;

% correction to solar longitude in arcsecs
slong_corr = -17.20 * sin(lanr) - 1.32 * sin(2.0*mslr) -0.23 * ...
    sin(2.0*mllr) + 0.21 * sin(2.0*lanr);

% conversion from arcsecs to degress
slong_corr = slong_corr/3600;

% correction to the obliquity
obliq_corr = 9.20 * cos(lanr) + 0.57 * cos(2.0*mslr) +0.10 * ...
    cos(2.0*mllr) - 0.09 * cos(2.0*lanr);

obliq_corr = obliq_corr / 3600;

last_jd = jd;
last_obliq_corr = obliq_corr;
last_slong_corr = slong_corr;

end





function mll = AstAlg_mean_lunar_longitude(jd)
%
% mll = AstAlg_mean_lunar_longitude(jd)
%
% Mean lunar longitude
%
% Based on the c-function AstAlg_mean_lunar_longitude in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  mll mean lunar longitude (deg)
%
% IV 2016

persistent last_jd last_mll

% if the value was already calculated, simply return it
if ~isempty(last_jd)
    if jd==last_jd;
        if ~isempty(last_mll)
            mll = last_mll;
            return
        end
    end
end

% time difference from J2000 in centuries...
tau = ( jd - J2000() ) / 36525;

% lunar longitude
llong = 218.3165 + 481267.8813 * tau;

mll = mod(llong,360);

last_jd = jd;
last_mll = mll;

end


function dec = AstAlg_solar_declination(jd)
%
% dec = AstAlg_solar_declination(jd)
%
% Solar declination
%
% Based on the c-function AstAlg_solar_declination in astalglib, see
% LICENSE-AstAlg.txt
%
% INPUT:
%  jd julian day
%
% OUTPUT:
%  dec solar declination (deg)
%
% IV 2016

persistent last_jd last_dec DTOR
if isempty(DTOR)
    DTOR = pi/180;
end

% if the value was already calculated, simply return it
if ~isempty(last_jd)
    if jd==last_jd;
        if ~isempty(last_dec)
            dec = last_dec;
            return
        end
    end
end

aob = AstAlg_apparent_obliquity(jd);
asl = AstAlg_apparent_solar_longitude(jd);

sindec = sin( DTOR * aob ) * sin( DTOR * asl );

dec = asin(sindec)/DTOR;

last_jd = jd;
last_dec = dec;

end