% Author: Tobias Gerken, tug15@psu.edu
% Version: 1.0 (2019_04_01)
% This code is used to generate ManeuverFlag files for ACT-America Flight
% Campaign using ACT-America Campaing Data:
%    ACT-America: L3 Merged In Situ Atmospheric Trace Gases and Flask
%    Data, Eastern USA", https://doi.org/10.3334/ORNLDAAC/1593

% ManeuverFlag definition: 0 = on ground; 
%                          1 = take off; 
%                          2 = in-line ascent; 
%                          3 = in-line descent; 
%                          4 = spiral up; 
%                          5 = spiral down; 
%                          6 = constant altitude legs; 
%                          7 = landing


% Needed: 
%   - ACT-America data netCDF files
%    Data product: "ACT-America: L3 Merged In Situ Atmospheric Trace Gases and Flask
%    Data, Eastern USA", https://doi.org/10.3334/ORNLDAAC/1593 
%    This algorithm uses the following contained variables (see ORNL documentation for details): 
%       'UTC', 'LATITUDE','LONGITUDE', 'GPS_ALT', 'Altitude-AGL_GoogleMaps', 'HDG'

% Output:
%   A CSV file per day containing the following data columns:
%   'UTC', 'LATITUDE', 'LONGITUDE', 'GPS_ALT'  - Values provided by ACT netCDF Data
%   'ManeuverFlag_Old' - Old ManeuverFlag variable as present in ACT netCDF Data
%   'ManeuverFlag_New' - New generated ManeuverFlag
%   'QC_MF' - A first quality control variable on Maneuver Flag; 
%           QC_MF = 0; probably low quality MF assignment or manucal
%                      adjustment of flag
%           QC_MF = 1; not data issue recorded (note this does not mean good data quality) 

% References/ Data Access:
%   Davis, K.J., M.D. Obland, B. Lin, T. Lauvaux, C. O'dell, B. Meadows, 
%   E.V. Browell, J.H. Crawford, J.P. DiGangi, C. Sweeney, M.J. McGill, 
%   J. Dobler, J.D. Barrick, and A.R. Nehrir. 2018. 
%   ACT-America: L3 Merged In Situ Atmospheric Trace Gases and Flask Data, 
%   Eastern USA. ORNL DAAC, Oak Ridge, Tennessee, USA. 
%   https://doi.org/10.3334/ORNLDAAC/1593
%  
% This code/script is provided for documentation purposes only and not endorsed by the ACT Project Team or NASA. 
% 
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
% INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH 
% THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Clear Memory
clear all; close all 
% From Matlab 2017a onwards legends autoupdate; prevent legends from
% updating on the fly. 
if ~verLessThan('matlab', '9.2')    
    set(0,'DefaultLegendAutoUpdate','off') 
end

%% Path defintion block
% These should be modified by user. 
basepath = './' ; % set base of path; all paths in script are relative to this base
codepath = [basepath 'ACTAmerica_ManeuverFlags/'] ; % add path to location of this code 
addpath(genpath(codepath)) ; 
% %%
 if strcmp(getenv('computername'),'DESKTOP-45CVB98')
     addpath(genpath('C:\Users\tobia\OneDrive - The Pennsylvania State University\Projects\ACT-America\Code\GitCode'));
     basepath = 'C:\Users\tobia\';
 elseif strcmp(getenv('computername'),'DESKTOP-A2GKIRA')
     addpath(genpath('D:\OneDrive - The Pennsylvania State University\Projects\ACT-America\Code\GitCode'));
     basepath = 'D:\';
 else
     addpath(genpath('C:\Users\tug15\OneDrive - The Pennsylvania State University\Projects\ACT-America\Code\GitCode'));
      basepath = 'C:\Users\tug15\';  
 end

codepath = [basepath 'OneDrive - The Pennsylvania State University\Projects\ACT-America\Code\ACT_ManeuverFlags\'] ; % add path to location of this code  
 
%% Settings
flag_debug = false; % Creates additional plots for debugging purposes 
set(0,'DefaultFigureVisible','on'); % create figures, but do not display if set to 'off'

DataPath.NC =  [basepath 'OneDrive - The Pennsylvania State University\Projects\ACT-America\Data\ACT_merge_share\current_version\'];
%DataPath.NC =  [codepath '\In\'] ; % Or any other path to netCDF files 
Dir.Out     =  [codepath 'Out\'] ;
mkdir(Dir.Out)

% ManeuverFlag definition: 0 = on ground; 1 = take off; 2 = in-line ascent; 3 = in-line descent; 4 = spiral up; 5 = spiral down; 6 = constant altitude legs; 7 = landing
FlagIsProf = [1 2 3 4 5 7];

% Code will run between Start- and EndDate 
StartDate = '2017-02-17';
EndDate   = '2018-08-31';

% Algorithm to calculate ManeuverFlag

% Defined thresholds
% Note: thresholds were set to achieve good agreement between expected
% Maneuver Flags and automatically assined Maneuver Flags. They don't neccesarily have
% an intrinsic physical meaning. 

threshold_ground = 25    ; % plane is considered on ground if Altitude agl is below (unit: m) 
threshold_landing = 375  ; % constraint on landing sequence
threshold_leg    = 300/5 ; % number of datapoints needed for a leg to be valid (5 minutes of flight time at present)
threshold_climb  = 5   ;   % minimum change in altitude needed for profile (unit: m), corresponds to 1 m/s at timsetp 5s
threshold_prof   = 500 ;   % minimum altitude difference (unit: m) needed for valid profile
threshold_high   = 925;    % level legs with mean(z_AGL) > threshold (unit: m) are considered high profiles 
                           % (z_MSL is governing level leg)
method_smooth ='movmedian';% select smoothing method to prevent very short segments
threshold_smooth = 20;     % select window for smoothing
threshold_lowprof = 750 ;  % threshold to calculate mean and standard deviation for low profiles  (unit: m)
noStdDev =  1.5;           % number of standard deviations (unit m) z_AGL needs to exceed for profile to start.
threshold_prominence = 200 ; % threshold for detection of local minima and maxima for profile detection

spiral_sectors    =12; % To test for a spiral profile a 360 degree rotation is divided into a number of sectors
spiral_threshold  =10; % Threshold for spiral detection; 
                       % E.g: For spiral_sectors = 12 and sprial_threshold = 10
                       % the plane heading must rotate through 300 degrees for
                       % profile to be considered spiral

% Quality control thresholds;                         
threshold_noisy  = 150; % Threshold for 'noisy' level leg determintation; corresponds to standard deviation of z_AGL 
                        % (unit m) 
threshold_buffer = 20 ; % Buffer (# datapoints) to prevent early cutoff of profiles
threshold_grad_2 = threshold_climb/2; % minimum change in altitude needed per timestep for noisy profiles 

%% Algorithm starts here; No modifications needed below

Planes = {'b200', 'c130'}; % Defines the Plane Names 
SDay = datenum(StartDate,'yyyy-mm-dd'); % Convert strings to iterable 
EDay = datenum(EndDate,'yyyy-mm-dd');

% Main loop over days and planes to assign ManeuverFlags
for day = SDay:EDay
fc(1:2)=1; 
for plane = Planes 
daystr = datestr(day,'yyyymmdd') % generate screen output for each day

if strcmp(plane,'b200')    
    % load netCDF Data for Flights  
    Vars =  {'UTC', 'LATITUDE','LONGITUDE', 'GPS_ALT','THETA', 'Altitude-AGL_GoogleMaps', 'HDG', ...
             'TEMPERATURE', 'TotalAirTemp',  'ManeuverFlag'};    
    f = dir([DataPath.NC 'ACTAMERICA-mrg05-b200_merge_' daystr '*.nc']);
    if ~isempty(f)
        info = ncinfo([DataPath.NC f.name]);
        NC = readNC_ORNL([DataPath.NC f.name], Vars);
    end
elseif strcmp(plane,'c130')
    Vars =  {'UTC', 'LATITUDE','LONGITUDE', 'GPS_ALT','THETA', 'Altitude-AGL_GoogleMaps', 'HDG'...
             'TEMPERATURE', 'TotalAirTemp', 'ManeuverFlag'};

    f = dir([DataPath.NC 'ACTAMERICA-mrg05-c130_merge_' daystr '*.nc']);
    if ~isempty(f)
        info = ncinfo([DataPath.NC f.name]);
        NC = readNC_ORNL([DataPath.NC f.name], Vars);
    end
end

if exist('NC','var') % skip days that don't have data 
% create temporary variables     
z_AGL = [NC.Altitude_AGL_GoogleMaps];
z_MSL = [NC.GPS_ALT];
dz_MSL = [diff(z_MSL) 0] ;
dz_AGL = [diff(z_AGL) 0] ;
Heading = [NC.HDG];
UTC = [NC.UTC]/3600;
MF = zeros(size(z_MSL))-1; % initialize ManeuverFlag with missing --> -1
MF_QC = zeros(size(z_MSL))+1; % initialize QCFlag with 1 --> good quality

% Core of algorithm begins here: 

% 0) set MF where altitude information is missing
MF(isnan(z_MSL)) = NaN;

% 1) Determine High Level Legs, during which Heigth MSL is constant
LL = abs(dz_MSL)< threshold_climb & z_AGL>threshold_high; 
% Eliminate short legs. A movmedian of 20 should remove legs with less than
% 30 seconds duration
LL = round(smoothdata(LL, method_smooth,threshold_smooth ));
LL = round(smoothdata(LL, method_smooth ,threshold_smooth )); % do twice to account for noisy first guess
% Assign level leg flag:
MF(logical(LL)) = 6 ;

% 2) Determine Low Level Legs, during which Heigth AGL more or less constant 
% Assume first that altitude agl is within 1 std of mean altitude below 750m
ind =  z_AGL<threshold_lowprof;
threshold_u = nanmean(z_AGL(ind))+noStdDev*nanstd(z_AGL(ind));
threshold_l = nanmean(z_AGL(ind))-noStdDev*nanstd(z_AGL(ind));
LL = z_AGL<threshold_u & z_AGL>threshold_l & z_AGL<=1000 & MF == -1 ; 
% Eliminate short legs. A movmedian of 20 should remove legs with less than
% 30 seconds duration
LL = round(smoothdata(LL, method_smooth,threshold_smooth ));
LL = round(smoothdata(LL, method_smooth,threshold_smooth )); % do twice to account for noisy first guess
% Assign level leg flag:
MF(logical(LL)) = 6 ; 

% 3) Determined airplane on ground
% Use a simple threshold, if plane is below 25m AGL it is considered on
% ground
ind = z_AGL < threshold_ground ;
MF(ind) = 0;

% First Quality control interlude
% The above algorith detects portions of profiles which are below 1000 m
% remove legs that have a mean gradient and that are shorter than 
D = (MF == 6);
FlagChange =[0 find(diff(D)) length(D)] ;
for i = 1:length(FlagChange)-1
    if MF(FlagChange(i)+1)~=6
        continue
    end
    if abs(nanmean(dz_AGL(FlagChange(i)+1:FlagChange(i+1))))>threshold_climb
       % remove legs that have mean gradient  
       MF(FlagChange(i)+1:FlagChange(i+1))=-1;
    elseif numel(FlagChange(i)+1:FlagChange(i+1)) < threshold_leg
       % remove short legs
       MF(FlagChange(i)+1:FlagChange(i+1))=-1;
    end   
end

% 4) Profiles;
% Loop over segments with unassigned ManeuverFlag and assign profile type flags
% MF = [1,2,3,4,5,7]
D = (MF == -1);
FlagChange =[0 find(diff(D)) length(D)];
for i = 1:length(FlagChange)-1 
    if MF(FlagChange(i)+1)~=-1 
        continue % ignore already assigned flags
    end
    % first determine wether we have a single profile
    % or several profiles in a row
    A = z_MSL(FlagChange(i)+1:FlagChange(i+1));
    t = UTC(FlagChange(i)+1:FlagChange(i+1));

    % find global minimum and maximum
    [~, glob_min_loc] = nanmin(A) ; 
    [~, glob_max_loc] = nanmax(A) ;
    % find local minima and maxima 
    loc_min_loc = find(islocalmin(A , 'MinProminence',threshold_prominence));
    loc_max_loc = find(islocalmax(A , 'MinProminence',threshold_prominence));
    
    ChangePoints = unique([1 sort([loc_min_loc loc_max_loc]) length(A)]);  
    
    if flag_debug
    figure
    x = 1:length(A);
    plot(t,A,t(ChangePoints),A(ChangePoints),'r*',t(loc_max_loc),A(loc_max_loc),'b*')
    end

    for p = 1:length(ChangePoints)-1       
        % Upwards legs 
        if A(ChangePoints(p)) < A(ChangePoints(p+1)) 
            is = max(FlagChange(i)+ChangePoints(p  ),1)          ;
            ie  = min(FlagChange(i)+ChangePoints(p+1),length(MF)) ; 
            iss = max(is-1, 1) ;
            if p == 1 && (MF(iss) == 0 |  FlagChange(i) ==0)
                % 4a) Takeoff
                MF(is:ie)=1;       
            else
                % ensure that heading rotates by 360; 
                % Each sector has to be present in profile
                x = floor(mod(Heading(is:ie),360)/(360/spiral_sectors)) ;
                s = unique(x);
                % Test weather spiral rotation exceeds minimum required
                % angle
                if ~(sum(ismember(s,[0:spiral_sectors-1])) >= spiral_threshold)
                    % Spiral condition not met
                    % 4b) Straight up
                    MF(is:ie) =  2;
                else
                    % Spiral condition met 
                    % 4c) Spiral up
                    MF(is:ie) =  4;
                end
            end
        else
            is = max(FlagChange(i)+ChangePoints(p  ),1)          ;
            ie = min(FlagChange(i)+ChangePoints(p+1),length(MF)) ;  
            iee = min(ie+1, length(MF)) ;
            if p == length(ChangePoints)-1 && (MF(iee) == 0 | FlagChange(i+1) == length(MF))
                % 4d) Landing 
                MF(is:ie)=7 ;
            
            else              
                % ensure that heading rotates by 360; 
                % Each sector has to be present in profile
                x = floor(mod(Heading(is:ie),360)/(360/spiral_sectors)) ;
                s = unique(x);
                
                if ~(sum(ismember(s,[0:spiral_sectors-1]))>= spiral_threshold)
                % Spiral condition not met:
                % 4e) Straight down
                    MF(is:ie) =  3;
                else
                    % Spiral condition met
                    % 4f) Spiral down
                    MF(is:ie) =  5;
                end
            end
        end
    end
    if flag_debug
        close(gcf)
    end
end

% temporarily remove NaN to allow for easier processing
MF(isnan(MF)) = -1;

% do a quality control loop over profile leg 
FlagChange = find(diff(MF));
for p = 1:length(FlagChange)-1
    fl = unique(MF(FlagChange(p)+1:FlagChange(p+1)));
        
    if ismember(MF(FlagChange(p)+1),[2:5]) && (numel(fl) == 1) 
        
        t = UTC(FlagChange(p)+1:FlagChange(p+1));
        z = z_AGL(FlagChange(p)+1:FlagChange(p+1));
        
        is = max(1,FlagChange(p)+1);
        ie = min(length(MF),FlagChange(p+1));        
        
        % Create figure for debugging
        if flag_debug
           figure
           hold on; 
           plot(t, z, 'color','red')
        end
        
        % extend profiles to nearest min/max in ABL.        
        if ismember(fl,[2,4]) && z(1) <1000 
           % rising profiles: look at values before 
           is = max(FlagChange(p)+1-threshold_buffer,1) ;
           t_ext = UTC(is:FlagChange(p));
           z_ext = z_AGL(is:FlagChange(p));
           
           % find location of minimum gradient near profile  
           g=smoothdata(gradient(z_ext));
           g_min = find(abs(g)==min(abs(g)));
           
           if flag_debug
              yyaxis left
              hold on
              plot(t_ext, z_ext, 'color','blue')
              plot(t_ext(g_min), z_ext(g_min), 'r*')
              yyaxis right
              grid on
              plot(t_ext, g)
           end
           
           % extend profile to minimum gradient, while preserving
           % surrounding profiles 
           is = is+g_min;
           ie = min(FlagChange(p+1),length(MF));
           ind = is:FlagChange(p);
           ind(MF(ind)~=6) = [] ;
           if ~isempty(ind)
              is = min(ind);
           else
              is = FlagChange(p)+1;
           end
           MF(ind) = fl ;
            
        elseif  ismember(fl,[3,5]) && z(end) <1000 
           % descending profile, look at values after  
           ie = min(FlagChange(p+1)+threshold_buffer,length(MF)) ;
           t_ext = UTC(FlagChange(p+1)+1:ie);
           z_ext = z_AGL(FlagChange(p+1)+1:ie);
           
           % find location of minimum gradient near profile 
           g=smoothdata(gradient(z_ext));
           g_min = find(abs(g)==min(abs(g)));
           
           if flag_debug
              yyaxis left
              hold on
              plot(t_ext, z_ext, 'color','blue')
              plot(t,z, 'color','red','linestyle','-')
              plot(t_ext(g_min), z_ext(g_min), 'r*')
              yyaxis right
              grid on
              plot(t_ext, g)
           end
           
           % extend profile to minimum gradient, while preserving
           % surrounding profiles 
           ie = ie+g_min;
           is = max(FlagChange(p)+1, 1);
           ind = FlagChange(p+1)+1:ie;
           ind(MF(ind)~=6) = [] ;
           if ~isempty(ind)
              ie = max(ind);
           else
              ie = FlagChange(p+1);
           end
           MF(ind) = fl ;
        end
        
        % do a quality check on profiles, to ensure that profiles have 
        % a) minimum range and total rise; b) minimum gradient; 
        r = range(z_AGL(is:ie));
        grad = abs(z_AGL(is) - z_AGL(ie)) / (ie-is);
        r2 = abs(z_AGL(is)-z_AGL(ie));
        
        if flag_debug
            text(0.02,0.9,['range ' num2str(r)],'units','normalized')
            text(0.02,0.8,['range2 ' num2str(r2)],'units','normalized')
            text(0.02,0.7,['grad ' num2str(grad)],'units','normalized')
        end
        
        % This is pretty arbitrarily chosen for now, but seems to mostly worked.   
        if ~(r > threshold_prof) || grad < threshold_grad_2 || ~(r2 > threshold_prof-100)
            MF(is:ie)=   6;    % set to leg for now
            MF_QC(is:ie)=0; % also assume that quality is less than perfect  
        end 
        if flag_debug
            close(gcf);
        end     
    end
end

FlagChange = find(diff(MF));

for p = 1:length(FlagChange)-1   
    fl = unique(MF(FlagChange(p)+1:FlagChange(p+1))) ;
    if ismember(MF(FlagChange(p)+1),7)   
    % do something
        if unique(MF(FlagChange(p-1)+1:FlagChange(p))) == 6
            zavg = nanmean(z_AGL(FlagChange(p-1)+1:FlagChange(p)));
            zavg = min(threshold_landing, zavg);
            if zavg<=threshold_high
                ind = find(z_AGL(FlagChange(p-1)+1:FlagChange(p))>= zavg,1,'last');
                MF(FlagChange(p-1)+ind:FlagChange(p)) = 7; 
            end
        end
    end
end

% Find level legs that have instances where z_AGL exceeds 1000m 
FlagChange = [0 find(diff(MF)) length(MF)];
ct = 0;
for p = 2:length(FlagChange)-1
    flags = MF(FlagChange(p)+1:FlagChange(p+1));
    if numel(unique(flags)) == 1  && unique(flags) == 6     
        ind = FlagChange(p)+1:FlagChange(p+1);
        % z_AGL(ind)
        avg = nanmean(z_AGL(ind));
        mx  = nanmax(z_AGL(ind)) ;
        if mx > threshold_high && avg < threshold_high 
            t= UTC(ind);
            z= z_AGL(ind);
            z_sm = smoothdata(z);
            
            if flag_debug
                figure
                box on; grid on;
                hold on
                h_dbg(1) = plot(t, z)     ;          
                h_dbg(2) = plot(t, z_MSL(ind));
                h_dbg(3) = plot(t, z_sm,'color','r');
                plot(t, t*0+avg,'color','k')   ;
            end
            
            if mx-avg > threshold_prof

                % find location of max
                i_mx = find(z==mx);
                % find location to left, where smoothed intersects avg
                i_1 = find(z_sm(1:i_mx)<=avg,1,'last');
                % find location to right, where smoothed profile intersects avg
                i_2 = find(z_sm(i_mx+1:end)<=avg,1,'first')+ i_mx;
                % determine whether profile condition is met and assign
                
                if flag_debug                    
                    h_dbg(4) = plot(t([i_1 i_mx i_2]),z_sm([i_1 i_mx i_2]),'r*') ;
                end
                
                if ~isempty(i_1)
                    MF([i_1:i_mx]+FlagChange(p)) = 2; % rising 
                    % assume that these are of low quality
                    MF_QC ([i_1:i_mx]+FlagChange(p)) = 0;
                end
                if ~isempty(i_2)
                    MF([i_mx+1:i_2]+FlagChange(p)) = 3; % descending
                    % assume that these are of low quality
                    MF_QC([i_mx+1:i_2]+FlagChange(p)) =0;
                end
            end
        
            if flag_debug
                legend(h_dbg,{'AGL','MSL','AGL_{sm}','Profs'},'location','best')
                ct = ct+1;
                print(gcf, [Dir.Out 'Debug_MtnProfs' daystr '_' plane{1} '_' num2str(ct) '.png'],'-dpng','-r300')
                clear h_dbg
                close(gcf)
            end
       
        end
    end
    
end

% Do a quality control over level legs 
% If level legs are too short they get merged with surrouding legs 
FlagChange = [0 find(diff(MF)) length(MF)];
for p = 2:length(FlagChange)-1
    
    flags = MF(FlagChange(p)+1:FlagChange(p+1));
    if numel(unique(flags)) == 1  && unique(flags) == 6  
    if numel(flags)<threshold_leg
        disp(numel(flags))
        % split up leg 
        ind = FlagChange(p)+1:FlagChange(p+1);
        MF_QC(ind) = 0; % assume that these are low quality
        ic = floor(length(ind)/2);
        
        mf_replace_b =  MF(FlagChange(p-1)+1:FlagChange(p));
        mf_replace_b = unique(mf_replace_b(isfinite(mf_replace_b))) ;
        mf_replace_a =  MF(FlagChange(p+1)+1:FlagChange(p+2)); 
        mf_replace_a  = unique(mf_replace_a(isfinite(mf_replace_a))) ;
        
        if ismember(mf_replace_b,[1, 2, 4]) && ismember(mf_replace_a,[3, 5, 7])
            MF(ind(1:ic))     = mf_replace_b ;
            MF(ind(ic+1:end)) = mf_replace_a ; 
        elseif ismember(mf_replace_a,[1, 2, 4]) && ismember(mf_replace_b,[3, 5, 7])
            MF(ind(1:ic))     = mf_replace_b ;
            MF(ind(ic+1:end)) = mf_replace_a ; 
        elseif mf_replace_b == 1 && ismember(mf_replace_a,[1, 2, 4])
            MF(ind) = 1;
            MF(FlagChange(p+1)+1:FlagChange(p+2)) = 1;
        elseif mf_replace_a == 7 && ismember(mf_replace_a,[3, 5, 7])
            MF(ind) = 7;
            MF(FlagChange(p-1)+1:FlagChange(p)) = 7;
        elseif mf_replace_a ==  mf_replace_b
            % merge legs of similar type
             MF(ind) = mf_replace_a;
            % do not replace if surrounding legs are in same direction but
            % of different type. 
        end

    end
    end
end

% Assign quality control over level legs
% If level legs are too short they get merged with surrouding legs 
FlagChange = [0 find(diff(MF)) length(MF)];
for p = 2:length(FlagChange)-1
    ind = FlagChange(p)+1:FlagChange(p+1);
    flags = MF(ind);
    if numel(unique(flags)) == 1  && unique(flags) == 6  
        z_avg = nanmean(z_AGL(ind));
        z_std = nanstd(z_AGL(ind)) ;
        if z_avg < threshold_high 
            if z_std> threshold_noisy
                if flag_debug
                   figure
                   hold on
                   plot(UTC(ind), z_AGL(ind))
                   plot([UTC(ind(1)) UTC(ind(end))],[z_avg+z_std z_avg+z_std],'-k',[UTC(ind(1)) UTC(ind(end))],[z_avg-z_std z_avg-z_std],'-k') 
                   close(gcf)
                end
                MF_QC(FlagChange(p)+1:FlagChange(p+1)) = 0;
            end            
        end
    end
        
end


MF(MF==-1) = NaN; % reset missing MF to NaN 
MF_QC(isnan(MF)) = NaN; % set MF_QC to NaN in MF is missing  

%% 
% Manually adjust some profiles. 
% I am doing this here rather than outside of Matlab for transparency and
% reproducability. 
% These are instances where the algorithm does not produce desired result
% and thus needs to be adjusted
if strcmp(plane,'b200')
    if     day == datenum('2016-07-21','yyyy-mm-dd')
        MF(59:121)    = 1;
        MF_QC(59:121) = 0;        
    elseif day == datenum('2016-07-25','yyyy-mm-dd')
        MF(2173:2259)    = 2;
        MF_QC(2173:2259) = 0;
    elseif day == datenum('2016-08-16','yyyy-mm-dd')
        MF(1796:1852)    = 7;
        MF_QC(1773:1796) = 0;
    elseif day == datenum('2016-08-21','yyyy-mm-dd')
        MF(6227:end)    = 7;
        MF_QC(6227:end) = 0;        
    elseif day == datenum('2017-02-17','yyyy-mm-dd')
        MF(4491:4522)    = 7;
        MF_QC(4491:4522) = 0;    
    elseif day == datenum('2017-11-03','yyyy-mm-dd')
        MF(1319:1368)    = 3;
        MF(1463:1496)    = 2;
        MF(1497:1605)    = 6;
    elseif day == datenum('2018-05-12','yyyy-mm-dd')
        MF_QC(1001:1110) = 0;
    end
elseif strcmp(plane,'c130')
    if     day == datenum('2016-08-22','yyyy-mm-dd') 
         MF(2938:3115)    = 2;
    elseif day == datenum('2017-03-09','yyyy-mm-dd') 
         MF(3849:3921)    = 6;
         MF(3921:3956)    = 7;
    elseif day == datenum('2017-03-10','yyyy-mm-dd')
         MF(2477:2513)    = 2;
         MF_QC(2477:2513) = 0;
    elseif day == datenum('2017-10-03','yyyy-mm-dd')
         MF(776:1057)    = 3;
         MF_QC(776:1057) = 0;
    elseif day == datenum('2017-10-14','yyyy-mm-dd')
         MF(635:655)    = 6;
    elseif day == datenum('2018-05-20','yyyy-mm-dd')
         MF(1749:1824)    = 3;
    end
end

if ~exist('fig','var')
    fig = figure('units','inches','position',[0 0 11 8.5]);
elseif ~ishandle(fig)
    fig = figure('units','inches','position',[0 0 11 8.5]);
end


if strcmp(plane,'b200')
    plt_start = 1;
    title_str = 'B200';
elseif strcmp(plane,'c130')
    plt_start = 2;
    title_str = 'C130';
end

cvec = cbrewer('qual','Set1',8);
colormap(cvec)
Leg = {'Ground', 'TakeOff', 'LegUp','LegDwn','SpiralU', 'SpiralD','Level','TouchD'};
subplot(3,2,plt_start)
box on; hold on; grid on;
title([title_str ' - ' datestr(day,'yyyy-mmm-dd')])
x=[NC.UTC]/3600;
scatter([[NC.UTC]/3600 x(end)+[1:8]/3600]  ,[[NC.GPS_ALT] [1:8]*0-100] ,20, [MF [0:7]],  'o',  'filled') 
ylabel('Height MSL')
ylm = get(gca,'ylim');
xlm = get(gca,'xlim');
set(gca,'ylim', [0 ylm(2)]) 
QC = MF_QC;
QC(MF_QC==1) = -50;
QC(QC==0) = 0.95*ylm(2);
plot([[NC.UTC]/3600 ],QC ,'linestyle','none','marker','.', 'linewidth',2, 'color','k') ;

if plt_start ==2 
%scatter([NC.UTC]/3600,[NC.Altitude_AGL_GoogleMaps], [], MF, 'marker','+', 'linewidth',2) 
    for t = 1:length(cvec)
    text(1.01,1-t*0.09, Leg{t},'color',cvec(t,:), 'fontsize',12,'fontweight','bold', 'units','normalized')
    end
end    
subplot(3,2,plt_start+2)
box on; hold on; grid on;
%scatter([NC.UTC]/3600,[NC.Altitude_AGL_GoogleMaps],20, MF,  'o',  'filled') 
% Quick and dirty fix for color map to make sure that all colors are
% present thus avoiding autmatic setting of coloraxis in matlab
scatter([[NC.UTC]/3600 x(end)+[1:8]/3600]  ,[[NC.Altitude_AGL_GoogleMaps] [1:8]*0-1000] ,20, [MF [0:7]],  'o',  'filled') 
ylabel('Height AGL')
ylm = get(gca,'ylim');
set(gca,'ylim', [0 ylm(2)]) 
set(gca,'xlim', xlm) 
QC = MF_QC;
QC(MF_QC==1) = -50;
QC(QC==0) = 0.95*ylm(2);
plot([[NC.UTC]/3600 ],QC ,'linestyle','none','marker','.', 'linewidth',2, 'color','k') ;


subplot(3,2,plt_start+4)
box on; hold on; grid on;
h(1) = plot([NC.UTC]/3600,[NC.ManeuverFlag],'linestyle','none','marker','.', 'linewidth',2, 'color',[0 0 0.75]) ;
h(2) = plot([NC.UTC]/3600,MF+0.2,'linestyle','none','marker','.', 'linewidth',2, 'color',[0.75 0 0]) ;
legend(h,{'netCDF','myAlg'},'location','northoutside','box','off','orientation','horizontal','fontsize',10)
set(gca, 'ylim', [-0.5 7.5], 'ytick',[0:7],'yticklabel',Leg)
set(gca,'xlim', xlm) 
ylabel('ManeuverFlag')

% WRITE TO CSV 
cHeader = {'UTC', 'LATITUDE','LONGITUDE', 'GPS_ALT', 'ManeuverFlag_Old' ,'ManeuverFlag_New', 'QC_MF' };
file = [Dir.Out 'ManeuverFlag_' plane{1} '_' daystr '.csv'];
data = [ [NC.UTC]; [NC.LATITUDE]; [NC.LONGITUDE]; [NC.GPS_ALT]; [NC.ManeuverFlag]; MF; MF_QC]' ;
[~] = writeCSV_head(file, data, cHeader);

% remove variables before next loop starts  
clear data file 
end

clear NC UTC dz_* z_* MF MF_Q

end

if exist('fig','var')
if ishandle(fig)
    print(fig, [Dir.Out 'ALTvMF' daystr '.png'],'-dpng','-r300') 
    close(fig)
end 
end
end
