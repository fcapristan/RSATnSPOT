%   LaunchSiteData.m: Compute derives launch site data and non-dimensionalize 
%
%   Michael Colonno, ADL, Stanford University 
%   Last Modified : Francisco Capristan (Added Vernal Equinox data)
%   Started:        7/1/11
%   Last Updated:   7/1/11
%
%   Inputs:         mission.launch_site         user launch site data   
%
%   Outputs:        mission.launch_site         updated / non-dimensional data 
%
%   CHANGES: Ross Allen
%   - Changed initial velocity to be in xy-plane instead of yz-plane

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mission = LaunchSiteData(mission)

mission.launch_site.non_dim.r = 1 + mission.launch_site.h/mission.scale.R;  % non-dim
lat = pi*mission.launch_site.lat/180;                                       % rad
long = pi*mission.launch_site.long/180;                                     % rad





local_time = mission.launch_site.local_time; % local launch time [hours, minutes]
local_time_to_UT = mission.launch_site.local_time_to_UT; %hour shift to convert local time to UT [hour]
year = mission.launch_site.date.year;
month = mission.launch_site.date.month;
day = mission.launch_site.date.day;

days = getDaysFromRefVernalEquinox(month,day,year,mission);

timeSEC = getUTsec(local_time,local_time_to_UT); %covert the local time [hours minutes] to UT seconds
D = days + timeSEC/(3600*24); % converting to fractions of DAYS ...ref Fund of ASTRO. Bates
thetag0 = mission.VernalEquinox.thetag0;
revPerDay = mission.VernalEquinox.revPerDay;
thetag = thetag0 + revPerDay*2*pi*D;

mission.launch_site.thetag = thetag;


%% x0 y0 z0 in IJK frame (see Vernal Equinox)
[x0, y0, z0] = sph2cart(long+thetag,lat,mission.launch_site.non_dim.r);            % non-dim
mission.launch_site.non_dim.r0 = [x0, y0, z0];

%% Using Earth's rotation at Planet's surface as default

% in V0 in IJK frame (see Vernal Equinox)    
mission.launch_site.non_dim.V0 = ...
    (mission.planet.omega*mission.planet.R*cos(lat)/mission.scale.V).*...
    [-sin(long+thetag), cos(long+thetag),0];                                             % non-dim

%% For Vertical or Horizontal Launches

if strcmp(mission.options.LaunchType,'vertical')
    rmag = sqrt(x0^2+y0^2+z0^2);
    u0 = [x0,y0,z0]/rmag;
    mission.launch_site.launchDirection = u0;
else
    error('Specify launch type. Currently only vertical launch supported');
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%