classdef Core_Sky < handle
    % This class contains properties and methods to manage astronomical objects
    
    %--- * --. --- --. .--. ... * ---------------------------------------------
    %               ___ ___ ___
    %     __ _ ___ / __| _ | __|
    %    / _` / _ \ (_ |  _|__ \
    %    \__, \___/\___|_| |___/
    %    |___/                    v 0.5.1 beta 3
    %
    %--------------------------------------------------------------------------
    %  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
    %  Written by: Giulio Tagliaferro
    %  Contributors:     Andrea Gatti
    %  A list of all the historical goGPS contributors is in CREDITS.nfo
    %--------------------------------------------------------------------------
    %
    %   This program is free software: you can redistribute it and/or modify
    %   it under the terms of the GNU General Public License as published by
    %   the Free Software Foundation, either version 3 of the License, or
    %   (at your option) any later version.
    %
    %   This program is distributed in the hope that it will be useful,
    %   but WITHOUT ANY WARRANTY; without even the implied warranty of
    %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %   GNU General Public License for more details.
    %
    %   You should have received a copy of the GNU General Public License
    %   along with this program.  If not, see <http://www.gnu.org/licenses/>.
    %
    %--------------------------------------------------------------------------
    % 01100111 01101111 01000111 01010000 01010011
    %--------------------------------------------------------------------------
    properties (Constant)
        group_delays_flags = [ 'GC1C' ; 'GC1S' ; 'GC1L' ; 'GC1X' ; 'GC1P' ; 'GC1W' ; 'GC1Y' ; 'GC1M' ; 'GC2C' ; 'GC2D' ; 'GC2S' ; 'GC2L' ; 'GC2X' ; 'GC2P' ; 'GC2W' ; 'GC2Y' ; 'GC2M' ; 'GC5I' ; 'GC5Q' ; 'GC5X' ; ... % GPS codes
            'RC1C' ; 'RC1P' ; 'RC2C' ; 'RC2P' ; 'RC3I' ; 'RC3Q' ; 'RC3X' ; ...                                                                         % GLONASS code
            'EC1A' ; 'EC1B' ; 'EC1C' ; 'EC1X' ; 'EC1Z' ; 'EC5I' ; 'EC5Q' ; 'EC5X' ; 'EC7I' ; 'EC7Q' ; 'EC7X' ; 'EC8I' ; 'EC8Q' ; 'EC8X' ; 'EC6A'; 'EC6B'; 'EC6C'; 'EC6X'; 'EC6Z';...          % GALILEO codes
            'JC1C' ; 'JC1S' ; 'JC1L' ; 'JC1X' ; 'JC1Z' ; 'JC2S' ; 'JC2L' ; 'JC2X' ; 'JC2M' ; 'JC5I' ; 'JC5Q' ; 'JC5X' ; 'JC6S' ; 'JC6L' ; 'JC6X' ; ... % QZSS codes
            'CC2I' ; 'CC2Q' ; 'CC2X' ; 'CC7I' ; 'CC7Q' ; 'CC7X' ; 'CC6I' ; 'CC6Q' ; 'CC6X' ; ...                                                       % BeiDou codes
            'IC5A' ; 'IC5B' ; 'IC5C' ; 'IC5X' ; 'IC9A' ; 'IC9B' ; 'IC9C' ; 'IC9X' ; ...                                                                % IRNSS codes
            'SC1C' ; 'SC5I' ; 'SC5Q' ; 'SC5X' % SBAS
            ];
    end
    
    
    properties
        time_ref_coord         % GPS Times of ephemerides
        time_ref_clock         %
        
        coord                  % Ephemerides [times x num_sat x 3]
        coord_type             % 0: Center of Mass 1: Antenna Phase Center
        poly_type              % 0: Center of Mass 1: Antenna Phase Center
        clock                  % clocks of ephemerides [times x num_sat]
        
        coord_rate = 900;
        clock_rate = 900;
        
        iono                   % 16 iono parameters
        
        X_sun                  % coord of sun ephemerides ECEF at the same time of coord
        X_moon                 % coord of moon ephemerides ECEF at the same time of coord
        sun_pol_coeff          % coeff for polynoimial interpolation of tabulated sun positions
        moon_pol_coeff         % coeff for polynoimial interpolation of tabulated moon positions
        
        erp                    % Earth Rotation Parameters
        
        % ALL Rinex 3 code observations flags + first letter indicationg the constellation
        
        group_delays = zeros(32,82); % group delay of code measurements (meters) referenced to their constellation reference:
        phase_delays = zeros(32,82); % group delay of code measurements (meters) referenced to their constellation reference:
        %    GPS     -> Iono free linear combination C1P C2P
        %    GLONASS -> Iono free linear combination C1P C2P
        %    Galileo -> Iono free linear combination
        %    BedDou  -> Iono free linear combination
        %    QZS     -> Iono free linear combination
        %    IRNSS   -> Iono free linear combination
        %    SABS    -> Iono free linear combination
        group_delays_times           % 77x1 GPS_Time
        
        tracking_bias
        
        ant_pco1              % satellites antenna phase center offset for the first frequency (used in toCOM toAPC)
        
        avail                 % availability flag
        coord_pol_coeff_apc   % coefficient of the polynomial interpolation for coordinates of the approximate antenna phase center  [11, 3, num_sat, num_coeff_sets]
        coord_pol_coeff_com   % coefficient of the polynomial interpolation for coordinates of the center of mass [11, 3, num_sat, num_coeff_sets]
        
        wsb                   % widelane satellite biases (only cnes orbits)
        wsb_date              % widelane satellite biases time
    end
    
    properties (Access = private)
        cc                    % constellation collector handler
    end
    
    methods
        % Creator
        function this = Core_Sky(force_clean)
            % Core object creator
            this.cc = Core.getCurrentSettings().getConstellationCollector;
            this.ant_pco1 = zeros(1, this.cc.getNumSat(), 3);
            if nargin == 1 && force_clean
                this.clearOrbit();
            end
        end
    end      
    
    % =========================================================================
    %  METHODS
    % =========================================================================
    
    methods % Public Access
        
        function is_empty = isEmpty(this)
            % Return true if the core_sky object does not contains data
            % At the momenth check for coord only
            %
            % SYNTAX
            %   sky.isEmpty
            is_empty = isempty(this.coord);
        end
        
        function initSession(this, start_date, stop_date, cc, flag_eph_only)
            % Load and precompute all the celestial parameted needed in a session delimited by an interval of dates
            % SYNTAX:
            %    this.initSession(this, start_date, stop_time)
            
            %%% load Epehemerids
            if nargin == 2 
                stop_date = start_date.last();
                start_date = start_date.first();
            end
            if nargin <= 3 || isempty(cc)
                new_cc = Core.getState.getConstellationCollector;
            else
                new_cc = cc;
            end
            flag_new_cc = ~strcmp(unique(this.cc.system), unique(new_cc.system));
            if flag_new_cc
                this.clearOrbit;
            end
            this.cc = new_cc;
            
            if nargin <= 4 || isempty(flag_eph_only)
                flag_eph_only = false;
            end
            
            flag_coo_loaded = ~isempty(this.getFirstEpochCoord) && this.getFirstEpochCoord <= start_date && this.getLastEpochCoord >= stop_date;
            flag_time_loaded = ~isempty(this.getFirstEpochClock) && this.getFirstEpochClock <= start_date && this.getLastEpochClock >= stop_date;
            %flag_update_sky = start_date < this.getCoordTime.first || stop_date > this.getCoordTime.last || ...
            %    start_date < this.getClockTime.first || stop_date > this.getClockTime.last;
            if ~isempty(start_date) && (flag_new_cc || (~flag_coo_loaded || (~flag_time_loaded && ~flag_eph_only)))
                this.group_delays = zeros(this.cc.getNumSat(),82); % group delay of code measurements (meters) referenced to their constellation reference:
                this.phase_delays = zeros(this.cc.getNumSat(),82); % group delay of code measurements (meters) referenced to their constellation reference:
                
                eph_f_name   = Core.getState.getEphFileName(start_date, stop_date);
                if isempty(eph_f_name)
                    flag_init_nav_files = true;
                else
                    [~, file_name, ~] = fileparts(eph_f_name{1});
                    if isempty(file_name)
                        flag_init_nav_files = true;
                    else
                        flag_init_nav_files = false;
                    end
                end
                
                if flag_init_nav_files
                    state = Core.getState();
                    center_name = state.getCurCenter();
                    if iscell(center_name)
                        center_name = center_name{1};
                    end
                    fw = File_Wizard;
                    list_preferred = state.PREFERRED_EPH(fw.rm.getOrbitType(center_name));
                    state.setPreferredOrbit(list_preferred, center_name);
                    eph_f_name   = Core.getState.getEphFileName(start_date, stop_date);
%                    fw.conjureNavFiles(start_date, stop_date);
                end
                
                clock_f_name = Core.getState.getClkFileName(start_date, stop_date);
                if isempty(clock_f_name)
                    clock_is_present = false;
                else
                    clock_is_present = true;
                    for i = 1:length(clock_f_name)
                        clock_is_present = clock_is_present && (exist(clock_f_name{i}, 'file') == 2);
                    end
                end
                clock_in_eph = isempty(setdiff(eph_f_name, clock_f_name)) || ~clock_is_present; %%% condition to be tested in different cases
                if isempty(this.time_ref_coord) || start_date < this.time_ref_coord
                    this.clearOrbit();
                else
                    to_clear_date = start_date.getCopy();
                    to_clear_date.addSeconds(-86400); % keep only one day before the first epoch
                    this.clearOrbit(to_clear_date);
                end
                
                this.poly_type = this.coord_type;
                if  instr(lower(eph_f_name{1}), '.sp3') || instr(lower(eph_f_name{1}), '.eph') || instr(lower(eph_f_name{1}), '.pre') % assuming all files have the same extensions
                    this.toCOM(); % be sure to be in COM before adding new coordinates
                    this.clearPolyCoeff();
                    this.clearSunMoon();
                    Core.getLogger.addMarkedMessage('Importing ephemerides...');
                    for i = 1:length(eph_f_name)
                        gps_time = getFileStTime(eph_f_name{i});
                        end_time = this.getLastEpochCoord();
                        % if true %isempty(end_time) || isempty(gps_time) || (end_time - gps_time) > -1e-3
                            this.addSp3(eph_f_name{i}, clock_in_eph);
                            %this.coord = this.coord(1 : find(any(this.coord(:,:,1),2), 1, 'last'),:,:);
                            %this.clock = this.clock(1 : find(any(this.clock(:,:),2), 1, 'last'),:,:);
                        % end
                        this.coord_type = 0; % center of mass
                        this.poly_type = 0; % center of mass
                    end
                    
                    % Very simple fill of missing orbits, but if they are bad, data are rejected
                    % in normal conditions orbits have no gaps, but sometimes BNC are not logged
                    % this part of code is commented because this "simple" fill will generate
                    % bad lagrange coefficients, a good interpolation have to pass probably for
                    % Keplerian parameters to be smoother
                    % for s = 1 : size(this.coord, 2)
                    %     if any(isnan(serialize(this.coord(:,s,:))))
                    %         for c = 1 : 3
                    %             data = this.coord(:,s,c);                                
                    %             id_ok = find(~isnan(data));
                    %             if numel(id_ok) > 3
                    %                 id_ko = find(isnan(data) & flagExpand(~isnan(data), 10)); % fill for a maximum of 10 epochs
                    %                 this.coord(id_ko,s,c) = interp1(id_ok, this.coord(id_ok,s,c), id_ko, 'spline', nan);
                    %             end
                    %         end
                    %     end
                    % end
                else %% if not sp3 assume is a rinex navigational file
                    clock_in_eph = true;
                    this.toAPC(); % be sure to be in APC before adding new coordinates
                    this.clearPolyCoeff();
                    this.clearSunMoon();
                    Core.getLogger.addMarkedMessage('Importing broadcast ephemerides...');
                    this.importBrdcs(eph_f_name,start_date, stop_date, clock_in_eph);
                    this.coord_type = 1; % antenna phase center
                    this.poly_type = 1;  % antenna phase center
                end
                Core.getLogger.addMarkedMessage('Precomputing satellite polynomial coefficients...');
                this.computeSatPolyCoeff();
                
                if Core.getState.isIonoKlobuchar
                    f_name = Core.getState.getIonoFileName(start_date, stop_date);
                    central_time = start_date.getCopy;
                    central_time.addSeconds((stop_date-start_date) / 2);
                    this.importIono(f_name{1}, central_time);
                end
                
                if not(clock_in_eph) && ~flag_eph_only
                    Core.getLogger.addMarkedMessage('Importing satellite clock files...');
                    for i = 1:length(clock_f_name)
                        [~,name,ext] = fileparts(clock_f_name{i});
                        gps_time = getFileStTime(eph_f_name{i});
                        end_time = this.getLastEpochClock();
                        if isempty(end_time) || isempty(gps_time) ||  gps_time > end_time
                            this.addClk(clock_f_name{i});
                            this.clock = this.clock(1 : find(any(this.clock(:,:),2), 1, 'last'),:,:);
                        end
                    end
                    
                    % Very simple fill of missing clocks, but if they are bad, data are rejected
                    % in normal conditions clocks have no gaps, but sometimes BNC are not logged
                    % NOTE CLOCK MUST NOT BE CHANGED -> very dangerous problems may arise
                    % This part of code have been commented because it is very difficult to
                    % fill the clocks in a good way. If no clocks are present => the
                    % satellite should not be used, but if the clocks are corrected then there
                    % is the risk to introduce an unhealty satellite causing bad estimations
                    % 
                    % for s = 1 : size(this.clock, 2)
                    %     if any(isnan(zero2nan(serialize(this.clock(:,s)))))
                    %         data = zero2nan(this.clock(:,s));
                    %         id_ok = find(~isnan(data));
                    %         if numel(id_ok) > 3
                    %             id_ko = find(isnan(data));
                    %             this.clock(id_ko,s) = nan2zero(interp1(id_ok, data(id_ok), id_ko, 'linear', nan));
                    %         end
                    %     end
                    % end
                end
                
                if not(flag_eph_only)
                    % Interp clock
                    this.fillClockGaps(10, 'spline'); % try to save small interval of missing clocks
                    
                    % load erp
                    Core.getLogger.addMarkedMessage('Importing Earth Rotation Parameters');
                    this.importERP(Core.getState.getErpFileName(start_date, stop_date),start_date);
                    
                    % load biases
                    Core.getLogger.addMarkedMessage('Importing code biases');
                    this.importBiases(Core.getState.getBiasFileName(start_date, stop_date));
                end
            end
        end
        
        function clearOrbit(this, gps_date)
            % clear the object of the data older than gps_date
            % SYNTAX: this.clearOrbit(gps_date)
            if nargin > 1
                this.clearCoord(gps_date);
                this.clearClock(gps_date);
            else
                this.clearCoord();
                this.clearClock();
            end
        end
        
        function clearCoord(this, gps_date)
            % DESCRIPTION: clear coord data, if date is provided clear
            % only data before that date
            if nargin > 1 && ~isempty(this.time_ref_coord)
                if this.time_ref_coord < gps_date
                    n_ep = min(floor((gps_date - this.time_ref_coord)/this.coord_rate), size(this.coord,1));
                    this.coord(1:n_ep,:,:)=[];
                    this.time_ref_coord.addSeconds(n_ep*this.coord_rate);
                    this.coord_pol_coeff_apc = []; %!!! the coefficient have to been recomputed
                    this.coord_pol_coeff_com = []; %!!! the coefficient have to been recomputed
                    
                    % delete also sun e moon data
                    if not(isempty(this.X_sun))
                        this.X_sun(1:n_ep,:)=[];
                    end
                    if not(isempty(this.X_moon))
                        this.X_moon(1:n_ep,:)=[];
                    end
                    this.sun_pol_coeff = []; %!!! the coefficient have to been recomputed
                    this.moon_pol_coeff = []; %!!! the coefficient have to been recomputed
                end
            else
                this.coord = [];
                this.time_ref_coord = [];
                this.coord_pol_coeff_apc = [];
                this.coord_pol_coeff_com = [];
                
                this.X_sun = [];
                this.X_moon = [];
                this.sun_pol_coeff = [];
                this.moon_pol_coeff = [];
            end
        end
        
        function clearClock(this, gps_date)
            % DESCRIPTION: clear clock data , if date is provided clear
            % only data before that date
            if nargin > 1  && ~isempty(this.time_ref_clock)
                if this.time_ref_clock < gps_date
                    n_ep = min(floor((gps_date - this.time_ref_clock)/this.clock_rate), size(this.clock,1));
                    this.clock(1:n_ep,:)=[];
                    this.time_ref_clock.addSeconds(n_ep*this.clock_rate);
                end
            else
                this.clock=[];
                this.time_ref_clock = [];
                this.wsb = [];
                this.wsb_date = [];
            end
        end
        
        function clearSunMoon(this, gps_date)
            % DESCRIPTION: clear sun and moon data , if date is provided clear
            % only data before that date
            if nargin > 1
                if this.time_ref_coord > gps_date
                    n_ep = floor((gps_date - this.time_ref_coord)/this.coord_rate);
                    this.X_sun(1:n_ep,:)=[];
                    this.X_moon(1:n_ep,:)=[];
                    this.sun_pol_coeff = []; %!!! the coefficient have to been recomputed
                    this.moon_pol_coeff = []; %!!! the coefficient have to been recomputed
                end
            else
                this.X_sun = [];
                this.X_moon = [];
                this.sun_pol_coeff = [];
                this.moon_pol_coeff = [];
            end
        end
        
        function clearPolyCoeff(this)
            % DESCRIPTION : clear the precomupetd poly coefficent
            this.coord_pol_coeff_apc = [];
            this.coord_pol_coeff_com = [];
            this.sun_pol_coeff = [];
            this.moon_pol_coeff = [];
        end
        
        function orb_time = getCoordTime(this)
            % DESCRIPTION:
            % return the time of coordinates in GPS_Time (unix time)
            orb_time = this.time_ref_coord.getCopy();
            orb_time.toUnixTime();
            [r_u_t , r_u_t_f ] = orb_time.getUnixTime();
            
            dt = (this.coord_rate : this.coord_rate : (size(this.coord,1)-1)*this.coord_rate)';
            
            
            u_t = r_u_t + uint32(fix(dt));
            u_t_f =  r_u_t_f  + rem(dt,1);
            
            idx = u_t_f >= 1;
            
            u_t(idx) = u_t(idx) + 1;
            u_t_f(idx) = u_t_f(idx) - 1;
            
            idx = u_t_f < 0;
            
            u_t(idx) = u_t(idx) - 1;
            u_t_f(idx) = 1 + u_t_f(idx);
            
            orb_time.appendUnixTime(u_t , u_t_f);
        end
        
        function orb_time = getClockTime(this)
            % DESCRIPTION:
            % return the time of clock corrections in GPS_Time (unix time)
            try
                orb_time = this.time_ref_clock.getCopy();
            catch
                Core.getLogger.addError('Core_Sky does not contains satellite clocks!');
                orb_time = [];
                return
            end
            orb_time.toUnixTime();
            
            [r_u_t , r_u_t_f ] = orb_time.getUnixTime();
            
            dt = (this.clock_rate : this.clock_rate : (size(this.clock,1)-1)*this.clock_rate)';
            
            u_t = r_u_t + uint32(fix(dt));
            u_t_f =  r_u_t_f  + rem(dt,1);
            
            idx = u_t_f >= 1;
            
            u_t(idx) = u_t(idx) + 1;
            u_t_f(idx) = u_t_f(idx) - 1;
            
            idx = u_t_f < 0;
            
            u_t(idx) = u_t(idx) - 1;
            u_t_f(idx) = 1 + u_t_f(idx);
            
            orb_time.appendUnixTime(u_t , u_t_f);
            
        end
        
        function time = getFirstEpochClock(this)
            % return first epoch of clock
            if ~isempty(this.time_ref_clock)
                time = this.time_ref_clock.getCopy();
            else
                time = [];
            end
        end
        
        function time = getLastEpochClock(this)
            % return last epoch of clock
            if ~isempty(this.time_ref_clock)
                time = this.time_ref_clock.getCopy();
                time.addSeconds((size(this.clock, 1) - 1) * this.clock_rate);
            else
                time = [];
            end
        end
        
        function time = getFirstEpochCoord(this)
            % return first epoch of coord
            if ~isempty(this.time_ref_coord)
                time = this.time_ref_coord.getCopy();
            else
                time = [];
            end
        end
        
        function time = getLastEpochCoord(this)
            % return last epoch of coord
            if ~isempty(this.time_ref_coord)
                time = this.time_ref_coord.getCopy();
                time.addSeconds((size(this.coord, 1) - 1) * this.coord_rate);
            else
                time = [];
            end
        end
        
        function eclipsed = checkEclipseManouver(this, time)
            eclipsed = int8(zeros(time.length,size(this.coord,2)));
            
            
            XS = this.coordInterpolate(time);
            
            %satellite geocentric position
            XS_n = sqrt(sum(XS.^2,3,'omitnan'));
            
            XS = XS./repmat(XS_n,1,1,3);
            
            %sun geocentric position
            X_sun = this.sunMoonInterpolate(time, true);
            X_sun = rowNormalize(X_sun);
            
            %satellite-sun angle
            cosPhi = sum(XS.*repmat(permute(X_sun,[1 3 2]),1,size(XS,2),1),3);
            %threshold to detect noon/midnight maneuvers
            thr = 4.9*pi/180*ones(time.length,size(this.coord,2)); % if we do not know put a conservative value
            
            shadowCrossing = cosPhi < 0 & XS_n.*sqrt(1 - cosPhi.^2) < (GPS_SS.ELL_A + 200000); % 200 Km buffer to be on the safe side
            atx = Core.getAntennaManager();
            sat_type = atx.getSatAntennaType(this.cc.getAntennaId, time.getCentralTime);
            if this.cc.isGpsActive
                for i = 1 : size(sat_type, 1) % only gps implemented
                    if length(sat_type{i}) > 8
                        if (strcmp(sat_type{i}(1:9),'BLOCK IIA'))
                            thr(:,i) = 4.9*pi/180; % maximum yaw rate of 0.098 deg/sec (Kouba, 2009)
                        elseif (strcmp(sat_type{i}(1:9),'BLOCK IIR'))
                            thr(:,i) = 2.6*pi/180; % maximum yaw rate of 0.2 deg/sec (Kouba, 2009)
                            % shadowCrossing(:,i) = false;  % shadow crossing affects only BLOCK IIA satellites in gps
                        elseif (strcmp(sat_type{i}(1:9),'BLOCK IIF'))
                            thr(:,i) = 4.35*pi/180; % maximum yaw rate of 0.11 deg/sec (Dilssner, 2010)
                            % shadowCrossing(:,i) = false;  % shadow crossing affects only BLOCK IIA satellites in gps
                        end
                    end
                end
            end
            %noon/midnight maneuvers affect all satellites
            noonMidnightTurn = acos(abs(cosPhi)) < thr;
            eclipsed(shadowCrossing) = 1;
            eclipsed(noonMidnightTurn) = 3;
        end
        
        function [is_shadowed] = isShadowed(this, XYZ, time)
            % compute weather the selected point is shadowed or not
            % considering earth as an ellipsoid
            %
            % SYNATX:
            % [is_shadowed] = isShadowed(this, XYZ, time)
            
            X_sun = this.sunMoonInterpolate(time, true);
            % find the angle over the equator
            alpha = atan2(X_sun(:,3), sqrt(sum(X_sun(:,1:2).^2,2)));
            lon_sun = atan2(X_sun(:,2),X_sun(:,1));
            % find earth radius at the tangent point (2d case)
            cc = Constellation_Collector('G');
            a = cc.gps.ELL_A;
            e = cc.gps.ELL_E;
            e2 = e^2;
            [Xi,Yi,Zi] = geod2cart(pi/2-alpha,lon_sun -pi,zeros(size(alpha)));
            b = sqrt(Xi.^2 +Yi.^2 +Zi.^2);
            % project tht epoint to be tested on the ellipse plane
            Z_rot = rowNormalize(X_sun);
            Y_rot = [Xi./b, Yi./b, Zi./b];
            X_rot = cross(Y_rot,Z_rot);
            XYZ_ell_plane = [sum(X_rot.*XYZ,2) sum(Y_rot.*XYZ,2) sum(Z_rot.*XYZ,2)];
            % check weather it is inside
            is_shadowed = (XYZ_ell_plane(:,1)./a).^2 + (XYZ_ell_plane(:,2)./b).^2 < 1 & XYZ_ell_plane(:,3) < 0;
            
        end
        
        function [sun_el, sun_az] = sunElevation(this, XYZ, time)
            % compute weather the selected point is shadowed or not
            % considering earth as an ellipsoid
            %
            % SYNATX:
            % [is_shadowed] = isShadowed(this, XYZ, time)
            
            X_sun = this.sunMoonInterpolate(time, true);
            [sun_az, sun_el] = deal(nan(size(XYZ,1),1));
            for i = 1 : size(XYZ,1)
                [sun_az(i), sun_el(i)] = this.computeAzimuthElevationXS(X_sun(i,:), XYZ(i,:));
            end
        end
        
        
        function importEph(this, eph, t_st, t_end, step, clock)
            % SYNTAX:
            %   eph_tab.importEph(eph, t_st, t_end, sat, step)
            %
            % INPUT:
            %   eph         = ephemerids matrix
            %   t_st        = start_time
            %   t_end       = start_time
            %   sat         = available satellite indexes
            %
            % OUTPUT:
            %   XS      = satellite position at time in ECEF(time_rx) (X,Y,Z)
            %   VS      = satellite velocity at time in ECEF(time_tx) (X,Y,Z)
            %   dtS     = satellite clock error (vector)
            %
            % DESCRIPTION:
            
            if nargin < 5
                step = 900;
            end
            this.coord_rate = step;
            if clock
                this.clock_rate = step;
            end
            if nargin < 4 || t_end.isempty()
                t_end = t_st;
            end
            times = (t_st.getGpsTime -5*step) : step : (t_end.getGpsTime+5*step); %%% compute 5 step before and after the day to use for polynomila interpolation
            this.time_ref_coord = t_st.getCopy();
            this.time_ref_coord.toUnixTime();
            this.time_ref_coord.addSeconds(-5*step);
            if clock
                this.time_ref_clock = this.time_ref_coord.getCopy();
            end
            this.coord = zeros(length(times), this.cc.getNumSat,3 );
            this.clock = zeros ( length(times),this.cc.getNumSat);
            systems = unique(eph(31,:));
            for sys = systems
                sat = unique(eph(30,eph(31,:) == sys)); %% keep only satellite also present in eph
                i = 0;
                prg_idx = sat; %this.cc.getIndex(sys,sat); % get progressive index of given satellites
                t_dist_exced = false;
                for t = times
                    i = i + 1;
                    [this.coord(i,prg_idx,:), ~, clock_temp, t_d_e, bad_sat] = this.satellitePositions(t, sat, eph(:, eph(31,:) == sys)); %%%% loss of precision problem should be less tha 1 mm
                    if clock
                        this.clock(i,prg_idx) = clock_temp';
                    end
                    t_dist_exced = t_dist_exced || t_d_e;
                end
                if t_dist_exced
                    cc = Core.getState.getConstellationCollector();
                    [ss, prn] = cc.getSysPrn(bad_sat);
                    str = '';
                    for s = 1 : numel(ss)
                        str = sprintf('%s, %c%02d', str, ss(s), prn(s));
                    end
                    str = str(3:end);
                    Core.getLogger.addWarning(sprintf('Satellite position problem:\nOne of the time bonds (%s , %s)\nfor sat %s\nis too far from valid ephemerids \nPositions might be inaccurate ',t_st.toString(0),t_end.toString(0), str));
                end
            end
        end
        
        function importBrdcs(this, f_names, t_st, t_end, clock, step)
            if nargin < 6
                step = 900;
            end
            if nargin < 5
                clock = true;
            end
            if nargin < 4 || t_end.isempty()
                t_end = t_st;
            end
            if not(iscell(f_names))
                f_names = {f_names};
            end
            eph = [];
            for i = 1:length(f_names)
                [eph_temp, this.iono] = this.loadRinexNav(f_names{i},this.cc,0,0);
                eph = [eph eph_temp];
            end
            
            if not(isempty(eph))
                this.importEph(eph, t_st, t_end, step, clock);
                %%% add TGD delay parameter
                for const = unique(eph(31,:))
                    eph_const = eph(:,eph(31,:)==const);
                    for s = unique(eph_const(1,:))
                        eph_sat = eph_const(:, eph_const(1,:) == s);
                        GD = eph_sat(28,1); % TGD change only every 3 months
                        
                        switch char(const)
                            case 'G'
                                idx_c1w = this.getGroupDelayIdx('GC1W');
                                idx_c2w = this.getGroupDelayIdx('GC2W');
                                this.group_delays(s,idx_c1w) = -GD * Core_Utils.V_LIGHT;
                                f = this.cc.getGPS().F_VEC; % frequencies
                                this.group_delays(s,idx_c2w) = - f(1)^2 / f(2)^2 * GD * Core_Utils.V_LIGHT;
                            case 'R'
                                idx_c1p = this.getGroupDelayIdx('RC1P');
                                idx_c2p = this.getGroupDelayIdx('RC2P');
                                this.group_delays(s,idx_c1p) = -GD * Core_Utils.V_LIGHT;
                                f = this.cc.getGLONASS().F_VEC; % frequencies
                                this.group_delays(s,idx_c2p) = - f(1)^2 / f(2)^2 * GD * Core_Utils.V_LIGHT;
                            case 'E'
                                idx_c1p = this.getGroupDelayIdx('EC1B');
                                idx_c2p = this.getGroupDelayIdx('EC5I');
                                this.group_delays(s,idx_c1p) = -GD * Core_Utils.V_LIGHT;
                                f = this.cc.getGalileo().F_VEC; % frequencies
                                this.group_delays(s,idx_c2p) = - f(1)^2 / f(2)^2 * GD * Core_Utils.V_LIGHT;
                            case 'C'
                                idx_c1p = this.getGroupDelayIdx('CC2I');
                                idx_c2p = this.getGroupDelayIdx('CC7I');
                                this.group_delays(s,idx_c1p) = -GD * Core_Utils.V_LIGHT;
                                f = this.cc.getBeiDou().F_VEC; % frequencies
                                this.group_delays(s,idx_c2p) = - f(1)^2 / f(2)^2 * GD * Core_Utils.V_LIGHT;
                                
                        end
                    end
                end
            end
        end
        
        function [XS,VS,dt_s, t_dist_exced, bad_sat] =  satellitePositions(this, time, sat, eph)
            
            % SYNTAX:
            %   [XS, VS] = satellitePositions(time_rx, sat, eph);
            %
            % INPUT:
            %   time_rx     = reception time
            %   sat         = available satellite indexes
            %   eph         = ephemeris
            %
            % OUTPUT:
            %   XS      = satellite position at time in ECEF(time_rx) (X,Y,Z)
            %   VS      = satellite velocity at time in ECEF(time_tx) (X,Y,Z)
            %   dtS     = satellite clock error (vector)
            %
            % DESCRIPTION:
            nsat = length(sat);
            
            XS = zeros(nsat, 3);
            VS = zeros(nsat, 3);
            
            
            dt_s = zeros(nsat, 1);
            t_dist_exced = false;
            bad_sat = [];
            for i = 1 : nsat
                
                k = find_eph(eph, sat(i), time, 86400);
                if not(isempty(k))
                    %compute satellite position and velocity
                    [XS(i,:), VS(i,:)] = this.satelliteOrbits(time, eph(:,k), sat(i), []);
                    dt_s(i) = sat_clock_error_correction(time, eph(:,k));
                    dt_s(i) = sat_clock_error_correction(time - dt_s(i), eph(:,k));
                else
                    t_dist_exced = true;
                    bad_sat = [bad_sat; sat(i)];
                end
                
            end
        end
        
        function addSp3(this, filename_SP3, clock_flag)
            % SYNTAX:
            %   this.addSp3(filename_SP3, clock_flag)
            %
            % INPUT:
            %   filename_SP3 = name of sp3 file
            %   clock_flag   = load also clock? (optional, default = true)
            %
            % DESCRIPTION:
            % add satellite and clock postiion contained in the sp3 file to
            % the object if values are contiguos with the ones already in
            % the object add them, otherwise clear the object and add them
            % data that are alrady present are going to be overwritten
            
            if isempty(this.coord)
                empty_file = true;
            else
                empty_file = false;
            end
            if nargin <3
                clock_flag = true;
            end
            
            % SP3 file
            f_sp3 = fopen(filename_SP3,'rt');
            
            if (f_sp3 == -1)
                Core.getLogger.addWarning(sprintf('No ephemerides have been found at %s the file cannot be opened', filename_SP3));
            else
                fnp = File_Name_Processor;
                Core.getLogger.addMessage(Core.getLogger.indent(sprintf('Opening file %s for reading', fnp.getFileName(filename_SP3))));
                
                txt = fread(f_sp3,'*char')';
                if isempty(txt)
                    Core.getLogger.addWarning(sprintf('No ephemerides have been found at "%s" the file is empty', filename_SP3));
                else
                    fclose(f_sp3);
                    if txt(1) == ' '
                        % the first line is corrupted and have spaces at
                        % the beginning of the line
                        id_cut = find(txt ~= ' ', 1, 'first');
                        if ~isempty(id_cut)
                            txt = txt(id_cut : end);
                        end
                    end
                    version = txt(2);
                    
                    % get new line separators
                    nl = regexp(txt, '\n')';
                    if nl(end) <  numel(txt)
                        nl = [nl; numel(txt)];
                    end
                    lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
                    lim = [lim lim(:,2) - lim(:,1)];
                    if lim(end,3) < 3
                        lim(end,:) = [];
                    end
                    % get end pf header
                    % coord  rate
                    coord_rate = cell2mat(textscan(txt(repmat(lim(2,1),1,11) + (26:36)),'%f'));
                    % n epochs
                    n_epochs = cell2mat(textscan(txt(repmat(lim(1,1),1,7) + (32:38)),'%f'));
                    % find first epoch
                    string_time = txt(repmat(lim(1,1),1,28) + (3:30));
                    
                    % import it as a GPS_Time obj
                    sp3_first_ep = GPS_Time(string_time, [], true);
                    if this.coord_rate ~= coord_rate
                        if empty_file
                            this.coord_rate = coord_rate;
                        else
                            Core.getLogger.addWarning(['Coord rate not match: ' num2str(coord_rate)]);
                            return
                        end
                    end
                    if clock_flag
                        this.clock_rate = coord_rate;
                    end
                    % checking overlapping and same correct syncro
                    sp3_last_ep = sp3_first_ep.getCopy();
                    sp3_last_ep.addSeconds(coord_rate * n_epochs);
                    
                    % initialize array size
                    if empty_file
                        this.time_ref_coord = sp3_first_ep.getCopy();
                        if clock_flag
                            this.time_ref_clock = sp3_first_ep.getCopy();
                        end
                        this.coord = zeros(n_epochs, this.cc.getNumSat(),3);
                        if clock_flag
                            this.clock = zeros(n_epochs, this.cc.getNumSat());
                        end
                    else
                        idx_first = (sp3_first_ep - this.time_ref_coord) / this.coord_rate;
                        idx_last = (sp3_last_ep - this.time_ref_coord) / this.coord_rate;
                        memb_idx = ismembertol([idx_first idx_last], -1 : (size(this.coord,1)+1) ); % check whether the extend of sp3 file intersect with the current data
                        c_n_sat = size(this.coord,2);
                        if memb_idx(1) == true && memb_idx(2) == false
                            n_new_epochs = idx_last - size(this.coord, 1);
                            if n_new_epochs > 0
                                this.coord = cat(1,this.coord,zeros(n_new_epochs,c_n_sat,3));
                                if clock_flag
                                    this.clock = cat(1,this.clock,zeros(n_new_epochs,c_n_sat));
                                end
                            end
                        elseif memb_idx(1) == false && memb_idx(2) == true
                            this.time_ref_coord = sp3_first_ep.getCopy();
                            if clock_flag
                                this.time_ref_clock = sp3_first_ep.getCopy();
                            end
                            n_new_epochs = -idx_first;
                            this.coord = cat(1,zeros(n_new_epochs,c_n_sat,3),this.coord);
                            if clock_flag
                                this.clock = cat(1,zeros(n_new_epochs,c_n_sat),this.clock);
                            end
                        elseif idx_first > 0 && (memb_idx(1) == false && memb_idx(2) == false)
                            % if there is a space between last coordinates and the new set
                            n_new_epochs = idx_last - size(this.coord, 1);
                            this.coord = cat(1,this.coord, zeros(n_new_epochs,c_n_sat,3));
                            if clock_flag
                                this.clock = cat(1,this.clock,zeros(n_new_epochs,c_n_sat));
                            end
                        end
                    end
                    %%%% read data
                    %%% read epochs
                    t_line = find(txt(lim(:,1)) == '*');
                    % find the length of the time string (it should of this format: yyyy mm dd HH MM SS.00000000)
                    % but sometimes it has 4 digits for seconds :-/ i.e. igs102664.sp3
                    % or has spaces at the end of the time line esa20341.sp3
                    time_len = length(strtrim(txt((lim(t_line(1),1)+1):lim(t_line(1),2)))) + 1;
                    string_time = txt(repmat(lim(t_line,1),1,time_len) + repmat(2 + (1:time_len), length(t_line), 1));
                    %string_time = [string_time repmat(' ',size(string_time,1),1)];
                    % import it as a GPS_Time obj
                    sp3_times = GPS_Time.fromString(string_time');
                    
                    % Get active constellations
                    % parse only these
                    sys_c_list = this.cc.getActiveSysChar;
                    % keep only valid lines
                    lim = lim(txt(lim(:,1)) == 'P' | txt(lim(:,1)) == '*', :);
                    
                    % Analyze data line and keep only the one with active satellites
                    d_line = find(txt(lim(:,1)) == 'P')';
                    sat_sys_c = txt(lim(d_line,1) + 1);
                    sat_sys_c(sat_sys_c == ' ') = 'G'; % standard "a" considers only GPS and no char G is written to specify the constellation
                    % keep only data lines with satellites that I want to load
                    id_ko = regexp(sat_sys_c, ['[^' sys_c_list ']']);
                    lim(d_line(id_ko), :) = [];
                    n_spe = diff(find([txt(lim(:, 1)) '*'] == '*')) - 1; % number of satellites per epoch
                    t_line = find(txt(lim(:, 1)) == '*');
                    d_line = find(txt(lim(:,1)) == 'P')';
                    
                    % Parse data considering 4 columns of valid observations with full size of 14chars (x4)
                    % All should be float values
                    data = reshape(sscanf(txt(bsxfun(@plus, repmat(lim(d_line,1) + 3, 1, 1 + 14*4), 1:(1 + 14*4)))', '%f'), 4, numel(d_line))';
                    
                    % character of sys_c for each line of the SP3
                    sys_c = txt(lim(d_line,1) + 1);
                    % list of epochs in this.coord to be written
                    c_ep_idx = round((sp3_times - this.time_ref_coord) / this.coord_rate) +1; %current epoch index
                    % prn of each line of data
                    sat_prn = uint8(txt(lim(d_line,1) + 2) - 48) * 10 + uint8(txt(lim(d_line,1) + 3) - 48);
                    
                    % id in this.coord containing the index to insert
                    id_ep = zeros(size(d_line,1),1);
                    id_ep(1) = 1;
                    id_next_ep = cumsum(n_spe(1 : end-1)) + 1; % id of the next epoch
                    for i = 1 : numel(id_next_ep)
                        id_ep(id_next_ep(i)) = id_ep(id_next_ep(i)) + 1;
                    end
                    id_ep = cumsum(id_ep);
                    
                    % go_id of each line of data
                    go_id = this.cc.getIndex(sys_c, double(sat_prn))';
                    
                    % remove unvalid PRNs
                    id_ko = isnan(go_id);
                    sat_prn(id_ko) = [];
                    data(id_ko,:) = [];
                    go_id(id_ko) = [];
                    d_line(id_ko) = [];
                    id_ep(id_ko) = [];
                    
                    % fill this.coord
                    for col = 1 : 3
                        id = c_ep_idx(id_ep) + (go_id - 1) * size(this.coord, 1) + (col - 1) * size(this.coord, 1) * size(this.coord, 2);
                        this.coord(id) = data(:, col) * 1e3;
                    end
                    
                    % fill this.clock
                    if clock_flag
                        data(:,4) = data(:,4) / 1e6;
                        data(data(:, 4) > 0.99, 4) = nan; % manual manage of nan (9999)
                        this.clock(c_ep_idx(id_ep) + (go_id - 1) * size(this.coord, 1)) = data(:, 4);
                    end
                end
            end
            this.coord = zero2nan(this.coord);  % <--- nan is slow for the computation of the polynomial coefficents
        end
                
        function fillClockGaps(this, max_gap, mode)
            % Fill clock gaps linearly interpolating neighbour clocks
            %
            % INPUT
            %   mode:   mode of interpolation
            %            - 'custom'    - linear custom made (default)
            %            - 'previous'  - Previous non-missing entry.
            %            - 'next'      - Next non-missing entry.
            %            - 'nearest'   - Nearest non-missing entry.
            %            - 'linear'    - Linear interpolation of non-missing entries.
            %            - 'spline'    - Piecewise cubic spline interpolation.
            %            - 'pchip'     - Shape-preserving piecewise cubic spline interpolation.
            %            - 'makima'    - modified Akima cubic interpolation.
            %
            % SYNTAX
            %   this.fillClockGaps(mode)
            %
            % SEE ALSO:
            %   fillmissing
            
            if nargin < 3
                mode = 'custom';
            end
            if ~exist('fillmissing', 'file')
                mode = 'custom';
            end
            
            switch mode
                case {'custom'}
                    for i = 1 : size(this.clock,2)
                        if not(sum(this.clock(:,i),1) == 0)
                            empty_clk_idx = this.clock(:,i) == 0 | isnan(this.clock(:,i));
                            clk_ok = ~empty_clk_idx;
                            n_ep = size(this.clock,1);
                            if sum(empty_clk_idx) < n_ep && sum(empty_clk_idx) > 0
                                this.clock(empty_clk_idx,i) = nan;
                                for hole = find(empty_clk_idx)'
                                    [idx_bf  ] = max((1 : hole)'   .* (this.clock(1 : hole ,i) ./this.clock(1 : hole ,i) ));
                                    [idx_aft ] = min((hole : n_ep)'.* (this.clock(hole : n_ep ,i) ./this.clock(hole : n_ep ,i)));
                                    if isnan(idx_bf)
                                        this.clock(hole,i) =  this.clock(idx_aft,i);
                                    elseif isnan(idx_aft)
                                        this.clock(hole,i) =  this.clock(idx_bf,i);
                                    else
                                        this.clock(hole,i) = ((idx_aft - hole) * this.clock(idx_bf,i) + (hole - idx_bf) * this.clock(idx_aft,i)) / (idx_aft - idx_bf);
                                    end
                                end
                            end
                            this.clock(~flagMergeArcs(clk_ok, max_gap), i) = nan;
                        end
                    end
                otherwise
                    for s = 1 : size(this.clock, 2)
                        % If there is a good observation and a missing value
                        if any(this.clock(:,s)) && any(isnan(zero2nan(this.clock(:,s))))
                            clk_ok = this.clock(:,s) ~= 0 & ~isnan(this.clock(:,s));                            
                            %this.clock(:,s) = fillmissing(zero2nan(this.clock(:,s)), mode, 'EndValues','nearest');
                            this.clock(:,s) = fillmissing(zero2nan(this.clock(:,s)), mode); % do not extrapolate
                            this.clock(~flagMergeArcs(clk_ok, max_gap), s) = nan;
                        end
                    end
            end
        end
        
        function addClk(this,filename_clk)
            % SYNTAX:
            %   eph_tab.addClk(filename_clk)
            %
            % INPUT:
            %   filename_clk = name of clk rinex file file (IMPORTANT:the method
            %   assume 1 day clock filen at 5s)
            %
            % DESCRIPTION:
            % add satellites  clock contained in the clk file to
            % the object if values are contiguos with the ones already in
            % the object add them, otherwise clear the object and add them
            % data that are alrady present are going to be overwritten
            f_clk = fopen(filename_clk,'rt');
            [~, fname, ~] = fileparts(filename_clk);
            if (f_clk == -1)
                Core.getLogger.addWarning(sprintf('No clk files have been found at %s', filename_clk));
            else
                fnp = File_Name_Processor;
                Core.getLogger.addMessage(Core.getLogger.indent(sprintf('Opening file %s for reading', fnp.getFileName(filename_clk))));
                t0 = tic;
                if isempty(this.clock)
                    empty_clk = true;
                else
                    empty_clk = false;
                end
                
                % read RINEX observation file
                try
                    txt = fread(f_clk,'*char')';
                catch ex
                    Core.getLogger.addWarning(sprintf('"%s" cannot be read correctly!', fnp.getFileName(filename_clk)));
                    Core_Utils.printEx(ex);
                end
                fclose(f_clk);
                
                % get new line separators
                nl = regexp(txt, '\n')';
                if nl(end) <  numel(txt)
                    nl = [nl; numel(txt)];
                end
                lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
                lim = [lim lim(:,2) - lim(:,1)];
                if lim(end,3) < 3
                    lim(end,:) = [];
                end
                
                % get end pf header
                eoh = strfind(txt,'END OF HEADER');
                eoh = find(lim(:,1) > eoh);
                eoh = eoh(1) - 1;
                if strcmp(fname(1:3),'grg') % if cnes orbit load wsb values
                    wl_line = txt(lim(1:eoh,1)) == 'W' & txt(lim(1:eoh,1)+1) == 'L'& txt(lim(1:eoh,1)+60) == 'C' & txt(lim(1:eoh,1)+61) == 'O' & txt(lim(1:eoh,1)+62) == 'M';
                    wsb_date = GPS_Time(cell2mat(textscan(txt(lim(find(wl_line,1,'first'),1) + [8:33]),'%f %f %f %f %f %f')));
                    wsb_prn = sscanf(txt(bsxfun(@plus, repmat(lim(wl_line, 1),1,3), 4:6))',' %f ');
                    wsb_value = sscanf(txt(bsxfun(@plus, repmat(lim(wl_line, 1),1,15), 39:53))','%f');
                    wsb = zeros(1,this.cc.getGPS.N_SAT);
                    wsb(wsb_prn) = wsb_value;
                    this.wsb = [this.wsb ;wsb];
                    this.wsb_date = [this.wsb_date ;wsb_date];
                end
                sats_line = find(txt(lim(eoh+1:end,1)) == 'A' & txt(lim(eoh+1:end,1)+1) == 'S') + eoh;
                % clk rate
                clk_rate = [];
                % find first epoch
                string_time = txt(repmat(lim(sats_line(1),1),1,27) + [8:34]);
                % convert the times into a 6 col time
                date = cell2mat(textscan(string_time,'%4f %2f %2f %2f %2f %10.7f'));
                % import it as a GPS_Time obj
                file_first_ep = GPS_Time(date, [], true);
                % find sampling rate
                ant_ids = this.cc.getAntennaId;
                ant_code_list = txt(bsxfun(@plus, repmat(lim(sats_line,1),1,3), 3:5));
                ant_id_list = Core_Utils.code3Char2Num(ant_code_list);
                for i = 1 : length(ant_ids)
                    ant_id = Core_Utils.code3Char2Num(ant_ids{i});
                    sat_line = sats_line(ant_id_list == ant_id);
                    if not(isempty(sat_line)) && length(sat_line) > 1
                        n_ep_sat = length(sat_line);
                        string_time = txt(repmat(lim(sat_line,1),1,27) + repmat(8:34, n_ep_sat, 1))';
                        % convert the times into a 6 col time
                        %date = cell2mat(textscan(string_time,'%4f %2f %2f %2f %2f %10.7f'));
                        % import it as a GPS_Time obj
                        %sat_time = GPS_Time(date, [], true);
                        sat_time = GPS_Time(string_time);
                        
                        % initilize matrix
                        if isempty(clk_rate)
                            clk_rate = median(diff(sat_time.getGpsTime()));
                            if not(empty_clk) & clk_rate ~= this.clock_rate
                                Core.getLogger.addWarning('Clock rate in file different from one in Core_Sky\n Discarding old data\n');
                                this.clearClock();
                                empty_clk = true;
                                this.clock_rate = clk_rate;
                            end
                            if empty_clk
                                this.clock_rate = clk_rate;
                                this.time_ref_clock = file_first_ep;
                                [ref_week, ref_sow] = this.time_ref_clock.getGpsWeek();
                                
                                this.clock = zeros(86400 / this.clock_rate, this.cc.getNumSat());
                            else
                                
                                c_ep_idx = round((file_first_ep - this.time_ref_clock) / this.clock_rate) +1; % epoch index
                                if c_ep_idx < 1
                                    this.clock = [zeros(abs(c_ep_idx)+1,size(this.clock,2)); this.clock];
                                    this.time_ref_clock = file_first_ep;
                                    [ref_week, ref_sow] =this.time_ref_clock.getGpsWeek();
                                end
                                c_ep_idx = round((file_first_ep - this.time_ref_clock) / this.clock_rate) +1; % epoch index
                                if c_ep_idx + 86400/this.clock_rate -1 > size(this.clock,1)
                                    this.clock = [this.clock; zeros( c_ep_idx + 86400/this.clock_rate -1 - size(this.clock,1) ,size(this.clock,2)); ];
                                end
                            end
                        end
                        c_ep_idx = round((sat_time - this.time_ref_clock) / this.clock_rate) +1; % epoch index
                        this.clock(c_ep_idx,i) = sscanf(txt(bsxfun(@plus, repmat(lim(sat_line, 1),1,21), 38:58))','%f');
                        if size(c_ep_idx, 1) ~= size(this.clock, 1) && Core.getState.isClockAlign()
                            % the clock is not empty and re-alignment have been requested
                            buf_size = 200;
                            % check left
                            id_left = max(1, c_ep_idx(1) - buf_size) : max(1, c_ep_idx(1) - 1);
                            id_left(this.clock(id_left, i) == 0) = [];
                            if numel(id_left) > 3
                                % realign on left
                                % compute id of the new clock set to be used for re-alignment
                                id_right = c_ep_idx(1) : c_ep_idx(find(c_ep_idx < (c_ep_idx(1) + buf_size), 1, 'last'));
                                id_right(this.clock(id_right, i) == 0) = [];
                
                                % DEBUG: inspect re-alignment
                                %figure(1000); clf; ax = axes;
                                %plot(ax, this.time_ref_clock.getMatlabTime + ((id_left(1) : id_right(end)) - 1) / 2880, Core_Utils.V_LIGHT * diff(this.clock(id_left(1):id_right(end)+1,i))); hold on;
                                %setTimeTicks
                                
                                if numel(id_right) > 3
                                    prediction_from_left = interp1(id_left, this.clock(id_left,i), c_ep_idx(1) + [-1 : 0], 'linear', 'extrap');
                                    prediction_from_right = interp1(id_right, this.clock(id_right,i), c_ep_idx(1) + [-1 : 0], 'linear', 'extrap');
                                    % figure(1000); clf; plot(id_left, this.clock(id_left,i), '.'); hold on; plot(id_right, this.clock(id_right,i), 'o');
                                    time_shift = mean(prediction_from_right - prediction_from_left);
                                    % Determine the treshold for jump detection
                                    % when the "jump" is lower than 5 times the std of the local change rate do not correct for it:
                                    thr = 1 * std([diff(zero2nan(this.clock(id_left(1) : id_left(end), i))); ...
                                                   diff(zero2nan(this.clock(id_right(1) : id_right(end),i)))], 'omitnan');
                                    if abs(time_shift) > thr 
                                        this.clock(c_ep_idx,i) = this.clock(c_ep_idx,i) - time_shift;
                                    %    plot([id_left id_right], this.clock([id_left id_right],i), 's')
                                    %else
                                    %    plot([id_left id_right], this.clock([id_left id_right],i), '^k')
                                    end
                                end
                                
                                % DEBUG: inspect re-alignment
                                %plot(ax, this.time_ref_clock.getMatlabTime + ((id_left(1) : id_right(end)) - 1) / 2880, Core_Utils.V_LIGHT * diff(this.clock(id_left(1):id_right(end)+1,i))); hold on;
                                %xlim(this.time_ref_clock.getMatlabTime + ([id_left(1) , id_right(end)] - 1) / 2880);
                                %ylabel('clock drift [m]');
                                %legend('clock error drift', 'clock error drift - re-aligned')
                                %fh = gcf; Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'light');                                
                            end
                            % check right
                            id_left = max(1, c_ep_idx(end) - buf_size + 1) : max(1, c_ep_idx(end));
                            id_left(this.clock(id_left, i) == 0) = [];
                            if numel(id_left) > 3
                                % realign on right
                                % compute id of the new clock set to be used for re-alignment
                                id_right = min(size(this.clock,1), (c_ep_idx(end) + 1)) : min(size(this.clock,1), (c_ep_idx(end) + buf_size));
                                id_right(this.clock(id_right, i) == 0) = [];
                                if numel(id_right) > 3
                                    prediction_from_left = interp1(id_left, this.clock(id_left,i), c_ep_idx(end) + [0 : 1], 'linear', 'extrap');
                                    prediction_from_right = interp1(id_right, this.clock(id_right,i), c_ep_idx(end) + [0 : 1], 'linear', 'extrap');
                                    time_shift = mean(prediction_from_right - prediction_from_left);
                                    % Determine the treshold for jump detection
                                    % when the "jump" is lower than 5 times the std of the local change rate do not correct for it:
                                    thr = 1 * std([diff(zero2nan(this.clock(id_left(1) : id_left(end), i))); ...
                                        diff(zero2nan(this.clock(id_right(1) : id_right(end),i)))], 'omitnan');
                                    if abs(time_shift) > thr
                                        this.clock((c_ep_idx(end) + 1): end,i) = this.clock((c_ep_idx(end) + 1) : end, i) - time_shift;
                                    end
                                end
                            end
                        end
                    end
                end
                Core.getLogger.addMessage(sprintf('Parsing completed in %.2f seconds', toc(t0)), 100);
                Core.getLogger.newLine(100);
                
                % Outlier detection
                %d_clock = Core_Utils.diffAndPred(zero2nan(this.clock));
                %d_clock = bsxfun(@minus, d_clock, median(d_clock, 'omitnan'));
                %d_clock = bsxfun(@minus, d_clock, median(d_clock, 2, 'omitnan'));
                %figure; plot(d_clock)
            end            
        end
        
        function importERP(this, f_name, time)
            this.erp = this.loadERP(f_name, time.getGpsTime());
        end
        
        function [erp, found] = loadERP(this, filename, time)
            % SYNTAX:
            %   [erp, found] = loadERP(filename, time);
            %
            % INPUT:
            %   filename = erp filename (including path) [string]
            %   time = GPS time to identify the time range of interest [vector]
            %
            % OUTPUT:
            %   erp = struct containing erp data
            %   found = flag to check if the required file was found
            %
            % DESCRIPTION:
            %   Tool for loading .erp files: Earth rotation parameters.
            
            fnp = File_Name_Processor();
            found = 0;
            erp = [];
            MJD = [];
            Xpole = [];
            Ypole = [];
            UT1_UTC = [];
            LOD = [];
            Xrt = [];
            Yrt = [];
            for f = 1 : length(filename)
                Core.getLogger.addMessage(Core.getLogger.indent(sprintf('Opening file %s for reading', fnp.getFileName(filename{f}))));
                fid = fopen(filename{f},'rt');
                
                if fid == -1
                    return
                end
                
                l=fgetl(fid);
                i=1;
                
                %check version
                if ~strcmp(l, 'version 2')
                    %wrong version
                    fclose(fid);
                    return
                end
                
                while isempty(strfind(l,'  MJD'));
                    if l==-1
                        fclose(fid);
                        return
                    end
                    l=fgetl(fid);
                    i=i+1;
                end
                i=i+1;
                fseek(fid, 0, 'bof');
                
                % [MJD,Xpole,Ypole,UT1_UTC,LOD,Xsig,Ysig,UTsig,LODsig,Nr,Nf,Nt,Xrt,Yrt,Xrtsig,Yrtsig] = textread(filename,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',' ','headerlines', i);
                
                ERP_data = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',i);
                
                if (isempty(ERP_data{1}))
                    fclose(fid);
                    return
                else
                    found = 1;
                end
                
                MJD = [MJD; ERP_data{1}]; %#ok<*AGROW>
                Xpole = [Xpole; ERP_data{2}];
                Ypole = [Ypole; ERP_data{3}];
                UT1_UTC = [UT1_UTC; ERP_data{4}];
                LOD = [LOD; ERP_data{5}];
                Xrt = [Xrt; ERP_data{13}];
                Yrt = [Yrt; ERP_data{14}];
                fclose(fid);
            end
            
            jd = MJD + 2400000.5;
            [gps_week, gps_sow, ~] = jd2gps(jd);
            [ERP_date] = gps2date(gps_week, gps_sow);
            [ERP_time] = weektow2time(gps_week, gps_sow ,'G');
            
            if ~any(ERP_time <= max(time) | ERP_time >= min(time))
                % no suitable epochs found in erp file
                erp = [];
                return
            end
            
            % assign erp values and compute rates (@ epoch of the first epoch of orbits
            % erp.t0 = min(time);
            
            %correct MJD with the length of day
            for i = 2 : length(ERP_time)
                ERP_time(i)=ERP_time(i)+(LOD(i-1)*0.1e-6);
            end
            
            erp.t = ERP_time;
            erp.Xpole = Xpole;
            erp.Ypole = Ypole;
            erp.Xrt = Xrt;
            erp.Yrt = Yrt;
            
            %coefficients of the IERS (2010) mean pole model
            t0 = 2000;
            cf_ante = [0   55.974   346.346; ...
                1   1.8243    1.7896; ...
                2  0.18413  -0.10729; ...
                3 0.007024 -0.000908];
            
            cf_post = [0   23.513   358.891; ...
                1   7.6141   -0.6287; ...
                2      0.0       0.0; ...
                3      0.0       0.0];
            
            idx_ante = find(ERP_date(:,1) <= 2010);
            idx_post = find(ERP_date(:,1)  > 2010);
            
            %computation of the IERS (2010) mean pole
            erp.meanXpole = zeros(size(erp.Xpole));
            erp.meanYpole = zeros(size(erp.Ypole));
            for d = 1 : 4
                if (~isempty(idx_ante))
                    erp.meanXpole(idx_ante) = erp.meanXpole(idx_ante) + (ERP_date(idx_ante,1) - t0).^cf_ante(d,1) .* cf_ante(d,2);
                    erp.meanYpole(idx_ante) = erp.meanYpole(idx_ante) + (ERP_date(idx_ante,1) - t0).^cf_ante(d,1) .* cf_ante(d,3);
                end
                
                if (~isempty(idx_post))
                    erp.meanXpole(idx_post) = erp.meanXpole(idx_post) + (ERP_date(idx_post,1) - t0).^cf_post(d,1) .* cf_post(d,2);
                    erp.meanYpole(idx_post) = erp.meanYpole(idx_post) + (ERP_date(idx_post,1) - t0).^cf_post(d,1) .* cf_post(d,3);
                end
            end
            
            erp.m1 =   erp.Xpole*1e-6 - erp.meanXpole*1e-3;
            erp.m2 = -(erp.Ypole*1e-6 - erp.meanYpole*1e-3);
        end
        
        function importBiases(this,fname)
            if not(isempty(fname))
                this.importSinexBias(fname{1});
            end
        end
        
        function importSinexBias(this,file_name)
            fid = fopen(file_name,'rt');
            if fid == -1
                Core.getLogger.addWarning(sprintf('Core_Sky: File %s not found', file_name));
                return
            end
            Core.getLogger.addMessage(Core.getLogger.indent(sprintf('Opening file %s for reading', file_name)));
            txt = fread(fid,'*char')';
            fclose(fid);
            
            % get new line separators
            nl = regexp(txt, '\n')';
            if nl(end) <  numel(txt)
                nl = [nl; numel(txt)];
            end
            lim = [[1; nl(1 : end - 1) + 1] (nl - 1)];
            lim = [lim lim(:,2) - lim(:,1)];
            if lim(end,3) < 3
                lim(end,:) = [];
            end
            
            % get end of header
            eoh = strfind(txt,'*BIAS SVN_ PRN ');
            eoh = find(lim(:,1) > eoh(1));
            
            eoh = eoh(1) - 1;
            head_line = txt(lim(eoh,1):lim(eoh,2));
            svn_idx = strfind(head_line,'PRN') - 1;
            c1_idx = strfind(head_line,'OBS1') -1 ;
            val_idx = strfind(head_line,'__ESTIMATED_VALUE____') - 1;
            std_idx = strfind(head_line,'_STD_DEV___') - 1;
            std_idx = std_idx(1);
            % removing header lines from lim
            lim(1:eoh, :) = [];
            
            % removing last two lines (check if it is a standard) from lim
            lim((end-1):end, :) = [];
            
            % removing non satellites related lines from lim
            sta_lin = txt(lim(:,1)+13) > 57 | txt(lim(:,1)+13) < 48; % Satellites have numeric PRNs
            lim(sta_lin,:) = [];
            
            % TODO -> remove bias of epoch different from the current one
            % find bias names presents
            fl = lim(:,1);
            
            tmp = [txt(fl+svn_idx)' txt(fl+svn_idx+1)' txt(fl+svn_idx+2)' txt(fl+c1_idx)' txt(fl+c1_idx+1)' txt(fl+c1_idx+2)'];
            idx = repmat(fl+val_idx,1,20) + repmat([0:19],length(fl),1);
            bias = sscanf(txt(idx)','%f');
            idx = repmat(fl+std_idx,1,11) + repmat([0:10],length(fl),1);
            bias_std = sscanf(txt(idx)','%f');
            % between C2C C2W the std are 0 -> unestimated
            % as a temporary solution substitute all the zero stds with the mean of all the read stds (excluding zeros)
            bad_ant_id = []; % list of missing antennas in the bias file
            bias_std(bias_std == 0) = mean(bias_std(bias_std ~= 0));
            ref_bias_name_old = '';
            bad_sat_str = '';
            for s = 1 : this.cc.getNumSat()
                sys = this.cc.system(s);
                prn = this.cc.prn(s);
                ant_id = this.cc.getAntennaId(s);
                sat_idx = this.prnName2Num(tmp(:,1:3)) == this.prnName2Num(ant_id);
                if sum(sat_idx) == 0
                    bad_ant_id = [bad_ant_id; ant_id];
                else
                    sat_bias_name = tmp(sat_idx,4:end);
                    sat_bias = bias(sat_idx);
                    try
                        sat_bias_std = bias_std(sat_idx);
                    catch ex
                        % no STD present
                    end
                    %check if there is the reference bias in the one
                    %provided by the external source
                    for b = 1 : size(sat_bias_name,1)
                        bias_col   = strLineMatch(this.group_delays_flags,[sys 'C' sat_bias_name(b,2:3)]);
                        if sat_bias_name(b,1) == 'C'
                            this.group_delays(s, bias_col) =  sat_bias(b) * Core_Utils.V_LIGHT * 1e-9;
                        else
                            this.phase_delays(s, bias_col) =  sat_bias(b) * Core_Utils.V_LIGHT * 1e-9;
                        end
                    end
                    
                end
            end
            if ~isempty(bad_sat_str)
                Core.getLogger.addWarning(sprintf('One or more biases are missing in "%s":\n%s\nthe bias will be eliminated only using iono-free combination', File_Name_Processor.getFileName(file_name), bad_sat_str));
            end
            if ~isempty(bad_ant_id)
                str = sprintf(', %c%c%c', bad_ant_id');
                if size(bad_ant_id, 1) > 1
                    Core.getLogger.addWarning(sprintf('Satellites %s not found in the bias file', str(3 : end)));
                else
                    Core.getLogger.addWarning(sprintf('Satellites %s not found in the bias file', str(3 : end)));
                end
            end
        end
        
        function idx = getGroupDelayIdx(this,flag)
            %DESCRIPTION: get the index of the gorup delay for the given
            %flag
            idx = find(sum(this.group_delays_flags == repmat(flag,size(this.group_delays_flags,1),1),2)==4);
        end
        
        function importIono(this, f_name, central_time)
            [~, this.iono, flag_return ] = this.loadRinexNav(f_name,this.cc, 1, Core.getCurrentSettings.getIonoModel, central_time);
            if (flag_return)
                return
            end
        end
        
        function [sx ,sy, sz] = getSatFixFrame(this, time, go_id)
            % SYNTAX:
            %   [i, j, k] = satellite_fixed_frame(time,X_sat);
            %
            % INPUT:
            %   time     = GPS_Time [nx1]
            %   X_sat    = postition of satellite [n_epoch x n-sat x 3]
            % OUTPUT:
            %   sx = unit vector that completes the right-handed system [n_epoch x n_sat x 3]
            %   sy = resulting unit vector of the cross product of k vector with the unit vector from the satellite to Sun [n_epoch x n_sat x 3]
            %   sz = unit vector pointing from the Satellite Mass Centre (MC) to the Earth's centre [n_epoch x n_sat x 3]
            %
            % DESCRIPTION:
            %   Computation of the unit vectors defining the satellite-fixed frame.
            
            
            t_sun = time;
            X_sun = this.sunMoonInterpolate(t_sun, true);
            if nargin > 2
                X_sat = this.coordInterpolate(time, go_id);
                X_sat = permute(X_sat,[1 3 2]);
            else
                X_sat = this.coordInterpolate(time);
            end
            n_sat = size(X_sat,2);
            sx = zeros(size(X_sat)); sy = sx; sz = sx;
            for idx = 1 : t_sun.length()
                x_sun = X_sun(idx,:);
                x_sat = X_sat(idx,:,:);
                e = permute(repmat(x_sun,1,1,n_sat),[1 3 2]) - x_sat ;
                e = e./repmat(normAlngDir(e,3),1,1,3); %sun direction
                k = -x_sat./repmat(normAlngDir(x_sat,3),1,1,3); %earth directions (z)
                j=cross(k,e); % perpendicular to bot earth and sun dorection (y)
                j= j ./ repmat(normAlngDir(j,3),1,1,3); % normalize, earth amd sun dorection are not perpendicular
                %                 j = [k(2).*e(3)-k(3).*e(2);
                %                     k(3).*e(1)-k(1).*e(3);
                %                     k(1).*e(2)-k(2).*e(1)];
                i=cross(j,k); %(x)
                %                 i = [j(2).*k(3)-j(3).*k(2);
                %                     j(3).*k(1)-j(1).*k(3);
                %                     j(1).*k(2)-j(2).*k(1)];
                sx(idx,:,:) = i ;
                sy(idx,:,:) = j ;
                sz(idx,:,:) = k ;
            end
            if n_sat == 1
                sx = squeeze(sx);
                sy = squeeze(sy);
                sz = squeeze(sz);
            end
            function nrm=normAlngDir(A,d)
                nrm=sqrt(sum(A.^2,d));
            end
        end
        
        function toCOM(this, force)
            %DESCRIPTION : convert coord to center of mass
            if ~isempty(this.coord)
                if isempty(this.coord_pol_coeff_apc)
                    this.computeSatPolyCoeff();
                end
                if isempty(this.coord_pol_coeff_com)
                    if this.coord_type == 1 % is APC
                        Core.getLogger.addMarkedMessage('Sat Ephemerids: switching to center of mass');
                        this.COMtoAPC(-1);
                        this.coord_type = 0;
                    end
                    this.computeSatPolyCoeff();
                end
                this.poly_type = 0;
            end
        end
        
        function toAPC(this, force)
            %DESCRIPTION : convert coord to center of antenna phase center
            if ~isempty(this.coord)
                if isempty(this.coord_pol_coeff_com)
                    this.computeSatPolyCoeff();
                end
                if isempty(this.coord_pol_coeff_apc)
                    if this.coord_type == 0
                        Core.getLogger.addMarkedMessage('Sat Ephemerids: switching to antenna phase center');
                        this.COMtoAPC(1);
                        this.coord_type = 1;
                    end
                    this.computeSatPolyCoeff();
                end
                this.poly_type = 1;
            end
        end
        
        function [wsb] = getWSB(this,time)
            [year_t, month_t, day_t ] = time.getCalEpoch();
            n_ep = size(this.wsb_date);
            n_ep = n_ep(1);
            wsb = 0;
            for i = 1:n_ep
                [year, month, day ] = this.wsb_date(i).getCalEpoch();
                if year == year_t && month == month_t && day == day_t
                    wsb = this.wsb(i,:);
                end
            end
        end
        
        function COMtoAPC(this, direction)
            if isempty(this.ant_pco1) || ~any(this.ant_pco1(:)) || size(this.ant_pco1, 2) ~= size(this.coord, 2)
                % load PCV
                Core.getLogger.addMarkedMessage('Loading antennas phase center offset');
                this.loadAntPCO();
            end
            [i, j, k] = this.getSatFixFrame(this.getCoordTime());
            sx = cat(3,i(:,:,1),j(:,:,1),k(:,:,1));
            sy = cat(3,i(:,:,2),j(:,:,2),k(:,:,2));
            sz = cat(3,i(:,:,3),j(:,:,3),k(:,:,3));
            % This is an approximation pco are different for various frequencies, in this way we move to the antenna center of L1!
            this.coord = this.coord + sign(direction)*cat(3, sum(repmat(this.ant_pco1, size(this.coord,1), 1, 1) .* sx , 3) ...
                , sum(repmat(this.ant_pco1,size(this.coord,1),1,1) .* sy , 3) ...
                , sum(repmat(this.ant_pco1,size(this.coord,1),1,1) .* sz , 3));
        end
                
        function pcv_delay = getPCV(this, ant, f_code, el, az)
            % Get the pcv correction for a given satellite antenna
            %
            % INPUT
            %   ant     Antenna object relative to the satellite
            %   f_code  frequency code for the PCV to extract
            %   el      elevation [deg]
            %   az      azimuth [deg]
            if this.poly_type == 0
                % if coordinates refers to center of mass apply also pco
                sat_pco = ant.getPCO(f_code) * 1e-3;
                neu_los = [cosd(az).*cosd(el) sind(az).*cosd(el) sind(el)];
                pco_delay = neu_los * sat_pco;
            else
                pco_delay = zeros(size(el));
            end
            pcv_delay = ant.getPCV(f_code, el, az) * 1e-3;            
            pcv_delay = pco_delay - pcv_delay;
        end

        function dts = getSatClock(this, gps_time, sat)
            % Get the clock of the satellite "sat" e.g. G21
            %
            % SYNTAX
            %   dts = this.getSatClock(gps_time <sat>)
            
            if nargin == 3
                go_id_list = this.cc.getIndex(sat);
                dts = this.clockInterpolate(gps_time, go_id_list);
            else
                dts = this.clockInterpolate(gps_time);
            end
        end
        
        function bias = getBias(this, go_id, o_code, gps_time)
            log = Core.getLogger();
            cc = Core.getConstellationCollector;
            if length(this.tracking_bias) >= go_id
                for b = 1 : length(this.tracking_bias{go_id})
                    if strcmpi(this.tracking_bias{go_id}{b}.o_code, o_code)
                        bias = this.tracking_bias{go_id}{b}.getBias(gps_time);
                        return
                    end
                end
                bias = zeros(gps_time.length,1);
                log.addWarning(sprintf('No bias correction present for obs %s sat %s',o_code,cc.getAntennaId(go_id)));
            else
                bias = nan(gps_time.length,1);
                log.addWarning(sprintf('No bias correction present for obs %s sat %s',o_code,cc.getAntennaId(go_id)));
            end
        end
        
        function [dts] = clockInterpolate(this, gps_time, sat_go_id)
            % SYNTAX:
            %   [dts] = clockInterpolate(time, sat);
            %
            % INPUT:
            %   time  = interpolation timespan GPS_Time
            %   sat   = satellite PRN
            %
            % OUTPUT:
            %   dts  = interpolated clock correction
            %
            % DESCRIPTION:
            %   SP3 (precise ephemeris) clock correction linear interpolation.
            if nargin < 3
                sat_go_id = this.cc.index;
            end
            
            dts = nan(gps_time.length, numel(sat_go_id));
            for s = 1 : numel(sat_go_id)
                sat = sat_go_id(s);
                interval = this.clock_rate;
                
                %find the SP3 epoch closest to the interpolation time
                %[~, p] = min(abs(SP3_time - time));
                % speed improvement of the above line
                % supposing SP3_time regularly sampled
                times = this.getClockTime();
                
                if isempty(times)
                    dts_tmp = 0;
                else
                    
                % find day change
                %date = times.get6ColDate;
                %day_change = find(diff(date(:,3)));
                
                p = max(1, min((round((gps_time - this.time_ref_clock) / interval) + 1)', times.length - 1));
                
                b =  (times.getEpoch(p) - gps_time)';
                
                SP3_c = zeros(gps_time.length,2);
                u = zeros(gps_time.length,1);
                
                % extract the SP3 clocks
                b_pos_idx = b > 0;
                p_pos = p(b_pos_idx);
                SP3_c(b_pos_idx,:) = cat(2, this.clock(max(p_pos - 1,1), sat), this.clock(p_pos, sat));
                u(b_pos_idx) = 1 - b(b_pos_idx)/interval;
                
                b_neg_idx = not(b_pos_idx);
                p_neg = p(b_neg_idx);
                SP3_c( b_neg_idx,:) = cat(2, this.clock(p_neg,sat), this.clock(p_neg+1,sat));
                u(b_neg_idx) = -b(b_neg_idx) / interval;
                
                dts_tmp = NaN * ones(size(SP3_c,1), size(SP3_c,2));
                idx = (sum(nan2zero(SP3_c) ~= 0,2) == 2 .* ~any(SP3_c >= 0.999,2)) > 0;
                
                %t = this.getClockTime.getRefTime(round(gps_time.first.getMatlabTime));
                data = this.clock(:, sat);
                %[~, ~, ~, dts_tmp] = splinerMat(t(not(isnan(data))), data(not(isnan(data))), 300, 1e-10, gps_time.getRefTime(round(gps_time.first.getMatlabTime)));
                dts_tmp = (1-u) .* SP3_c(:,1) + (u) .* SP3_c(:,2);
                % detect anomalous clock jumps ( > 5 cm) if the clack have
                % a rate of at least 5 minutes
                if this.getClockTime.getRate >= 300
                    id_ko = find(abs(diff(data) - median(diff(data))) * Core_Utils.V_LIGHT > 0.05);
                    for i = 1 : numel(id_ko)
                        id_ko_tmp = find((p == id_ko(i) & b_neg_idx) | ...
                            (p == id_ko(i) + 1 & b_pos_idx));
                        %lid_ok = flagExpand(lid_ko, 5);
                        if numel(id_ko_tmp) > 3
                            % This will force to split the system
                            dts_tmp(id_ko_tmp([2 end-1])) = NaN;
                        end
                    end
                end
                clear data;
                dts_tmp(not(idx)) = NaN;
                
                %             dt_S_SP3=NaN;
                %             if (sum(SP3_c~=0) == 2 && ~any(SP3_c >= 0.999))
                %
                %                 %linear interpolation (clock)
                %                 dt_S_SP3 = (1-u)*SP3_c(1) + u*SP3_c(2);
                %
                %                 %plot([0 1],SP3_c,'o',u,dt_S_SP3,'.')
                %                 %pause
                %             end
                end
                if numel(sat_go_id) == 1
                    dts = dts_tmp(:);
                else
                    dts(:,s) = dts_tmp(:);
                end
            end
            
        end
        
        function computeSatPolyCoeff(this, order, n_obs)
            % SYNTAX:
            %   this.computeSatPolyCoeff(degree, n_obs);
            %
            % INPUT:
            %
            % OUTPUT:
            %
            % DESCRIPTION: Precompute the coefficient of the N th poynomial for all the possible support sets
            if nargin == 1
                order = 10;
                n_obs = 11;
            end
            order = order + mod(order, 2); % order needs to be even
            n_obs = max(order + 1, n_obs); % minimum num of observations to use is == num of coefficients
            n_obs = n_obs + mod(n_obs + 1, 2); % n_obs needs to be odd
            
            n_coeff = order + 1;
            A = zeros(n_obs, n_coeff);
            A(:, 1) = ones(n_obs, 1);
            x = -((n_obs - 1) / 2) : ((n_obs - 1) / 2); % * this.coord_rate
            for i = 1 : order
                A(:, i + 1) = (x .^ i)';         % y = a + b*x + c*x^2
            end
            n_coeff_set = size(this.coord, 1) - n_obs + 1; %86400/this.coord_rate+1;
            n_sat = size(this.coord, 2);
            coord_pol_coeff = zeros(n_coeff, 3, n_sat, n_coeff_set);
            % There is an equivalent implementation but a lot faster => see below
            % Keep these lines as a reference because results are slightly different (numerical differences)
            %for s = 1 : n_sat
            %    for i = 1 : n_coeff_set
            %        for j = 1 : 3
            %            % this.coord_pol_coeff(: , j, s, i) = (A' * A) \ A' * squeeze(this.coord(i : i + n_obs - 1, s, j));
            %            coord_pol_coeff(: , j, s, i) = A \ squeeze(this.coord(i : i + n_obs - 1, s, j));
            %        end
            %    end
            %end
            % Equivalent implementation but a lot faster
            tmp_coord = permute(this.coord, [1 3 2]);
            c_size = size(tmp_coord);
            a_size = size(A);
            inv_A = inv(A);
            for i = 1 : n_coeff_set
                coord_pol_coeff(: , :, :, i) = reshape(inv_A * reshape(tmp_coord(i : i + n_obs - 1, :, :), a_size(1), c_size(2) * c_size(3)), a_size(1), c_size(2), c_size(3));
            end
            
            if this.coord_type == 1 % APC
                this.coord_pol_coeff_apc = coord_pol_coeff;
            else % COM
                this.coord_pol_coeff_com = coord_pol_coeff;
            end
        end
        
        function computeSat11PolyCoeff(this)
            % SYNTAX:
            %   this.computeSat11PolyCoeff();
            %
            % INPUT:
            %
            % OUTPUT:
            %
            % DESCRIPTION: Precompute the coefficient of the 10th poynomial for all the possible support sets
            n_pol = 10;
            n_coeff = n_pol + 1;
            A = zeros(n_coeff, n_coeff);
            A(:, 1) = ones(n_coeff, 1);
            x = -5 : 5; % *this.coord_rat
            for i = 1 : 10
                A(:, i + 1) = (x .^ i)';
            end
            n_coeff_set = size(this.coord, 1) - 10; %86400/this.coord_rate+1;
            n_sat = size(this.coord, 2);
            coord_pol_coeff = zeros(n_coeff, 3, n_sat, n_coeff_set);
            for s = 1 : n_sat
                for i = 1 : n_coeff_set
                    for j = 1 : 3
                        coord_pol_coeff(: , j, s, i) = A \ squeeze(this.coord(i : i + 10, s, j));
                    end
                end
            end
            
            if this.coord_type == 1 % APC
                this.coord_pol_coeff_apc = coord_pol_coeff;
            else % COM
                this.coord_pol_coeff_com = coord_pol_coeff;
            end
        end
        
        function computeSMPolyCoeff(this)
            % SYNTAX:
            %   this.computeSMPolyCoeff();
            %
            % INPUT:
            %
            % OUTPUT:
            %
            % DESCRIPTION: Precompute the coefficient of the 10th poynomial for all the possible support sets
            n_pol = 10;
            n_coeff = n_pol+1;
            A = zeros(n_coeff, n_coeff);
            A(:, 1) = ones(n_coeff, 1);
            x = -5 : 5; % *this.coord_rat
            for i = 1 : 10
                A(:,i+1) = (x.^i)';
            end
            n_coeff_set = size(this.X_sun, 1) - 10;%86400/this.coord_rate+1;
            this.sun_pol_coeff = zeros(n_coeff, 3, n_coeff_set);
            this.moon_pol_coeff = zeros(n_coeff, 3, n_coeff_set);
            for i = 1 : n_coeff_set
                for j = 1 : 3
                    this.sun_pol_coeff(:,j,i) = A \ squeeze(this.X_sun(i:i+10,j));
                    this.moon_pol_coeff(:,j,i) = A \ squeeze(this.X_moon(i:i+10,j));
                end
            end
        end
        
        function [coord_pol_coeff, flag_update] = getPolyCoeff(this)
            % Get polynomial coefficients of the orbits
            %
            % SYNTAX
            %   [coord_pol_coeff, flag_update] = this.getPolyCoeff()
            flag_update = false;
            
            if this.poly_type == 1 % APC
                if isempty(this.coord_pol_coeff_apc)
                    flag_update = true;
                    if this.coord_type == 0 % COM
                        this.toAPC
                    else
                        this.computeSatPolyCoeff();
                    end
                end
                coord_pol_coeff = this.coord_pol_coeff_apc;
            else % COM
                if isempty(this.coord_pol_coeff_com)
                    flag_update = true;
                    if this.coord_type == 1 % APC
                        this.toCOM
                    else
                        this.computeSatPolyCoeff();
                    end
                end
                coord_pol_coeff = this.coord_pol_coeff_com;
            end
        end
                
        function [az, el, sat_coo] = getAzimuthElevation(this, coo, gps_time, sat)
            % Get the azimuth and elevation of a satellite "sat" e.g. G21
            % 
            % INPUT
            %   gps_time    GPS_Time object
            %   sat         satellite (eg. 'R02') / go_id (e.g. 34)
            %
            % OUTPUT
            %   az, el      azimuth and elevation [deg]
            %   sat_coo     coordinates of the satellite (Coordinate object)
            %
            % SYNTAX
            %   [az, el, sat_coo] = getAzimuthElevation(this, coo, gps_time, sat)
            
            if nargin == 4
                sat_xyz = this.getSatCoord(gps_time, sat);
            else
                sat_xyz = this.getSatCoord(gps_time);
            end
            [az, el] = this.computeAzimuthElevationXS(sat_xyz, coo.getXYZ());
            
            if nargout == 3
                sat_coo = Coordinates.fromXYZ(sat_xyz);
            end
        end
        
        function [X_sat, V_sat] = getSatCoord(this, gps_time, sat)
            % Get the coordinate of the satellite "sat" e.g. G21
            %
            % INPUT
            %   gps_time    GPS_Time object
            %   sat         satellite (eg. 'R02') / go_id (e.g. 34)
            %
            % SYNTAX
            %   [X_sat, V_sat] = this.getSatCoord(gps_time, <sat>)
            
            if nargin == 3
                if isnumeric(sat)
                    go_id_list = sat;
                else
                    go_id_list = this.cc.getIndex(sat);
                end
                [X_sat, V_sat] = this.coordInterpolate(gps_time, go_id_list);
            else
                [X_sat, V_sat] = this.coordInterpolate(gps_time);
            end
        end        
        
        function coordFit(this, gps_time, coord, go_ids)
            % Fit coordinates of satellites and add them to current
            % estimates
            %
            % INPUT:
            %    gps_time = time of the coordinates
            %    coord = value of the coordinates
            %    go_ids = go id of the satellites
            %
            % SYNTAX:
            %   this.coordFit(gps_time, coord, go_ids)
            
            %number of seconds in a quarter of an hour
            interval = this.coord_rate;
            n_obs = gps_time.length;
            %find the SP3 epoch closest to the interpolation time
            %[~, p] = min(abs(SP3_time - time));
            % speed improvement of the above line
            % supposing SP3_time regularly sampled
            t_diff = gps_time.getRefTime(this.time_ref_coord.getMatlabTime);

            p = round((t_diff / interval)+eps(t_diff / interval)) + 1;
            
            b = (p-1)*interval - t_diff;
                        
            %Lagrange interpolation
            %degree of interpolation polynomial (Lagrange)

            u = 6 + (- b )/interval;    % using 6 since n = 10;
            A = Core_Sky.fastLI(u)';
            A_idx = repmat(p,1,11) + repmat((-5:5),n_obs,1);
            idx = A_idx > 0;
            rows = repmat((1:n_obs)',1,11);
            n_par = max(max(A_idx));
            A = sparse(rows(idx),A_idx(idx),A(idx),n_obs,n_par);
            idx_empty = sum(A~=0,1) == 0;
            A(:,idx_empty) = [];
            n_par = size(A,2);
            N = A'*A + eye(n_par)*0.001;
            n_coord = size(this.coord,1);
            for s = 1 : length(go_ids)
                for c = 1 : 3
                    B = A'*coord(:,s,c);
                    coord_cord = zeros(n_coord,1);
                    coord_cord(~idx_empty) = N \ B;
                    this.coord(:,go_ids(s),c) = this.coord(:,go_ids(s),c) + coord_cord(1:min(n_coord,length(coord_cord)));
                end
            end
            
        end
        
        
        function [X_sat, V_sat] = coordInterpolate(this, gps_time, go_id)
            % Interpolate coordinates of satellites
            %
            % INPUT:
            %    t = vector of times where to interpolate
            %    sat = satellite to be interpolated (optional)
            %
            % SYNTAX:
            %   [X_sat] = Eph_Tab.polInterpolate(t, sat)
            
           if isempty(this.time_ref_coord)
                Core.getLogger.addWarning('Core_Sky appears to be empty, goGPS is going to miesbehave\nTrying to load needed data')
                this.initSession(gps_time.first(), gps_time.last())
            end
            n_sat = size(this.coord, 2);
            if nargin <3
                sat_idx = ones(n_sat, 1) > 0;
            else
                sat_idx = go_id;
            end
            
            n_sat = length(sat_idx);
            poly = this.getPolyCoeff;
            poly_order = size(poly,1) - 1;
            n_border = ((size(this.coord, 1) - size(poly, 4)) / 2);
            
            % Find the polynomial id at the interpolation time
            %t_diff = round(t.getRefTime(this.time_ref_coord.getMatlabTime),7); % round to 7 digits, t - ref_time cannot hold more precision
            t_diff = gps_time.getRefTime(this.time_ref_coord.getMatlabTime);
            
            pid_floor = floor(t_diff / this.coord_rate) + 1 - n_border;
            % Ignore solution at the border of the polynomial
            pid_floor(pid_floor < 1) = 1;
            pid_floor(pid_floor > size(this.getPolyCoeff, 4)) = size(this.getPolyCoeff, 4);
            pid_ceil = ceil(t_diff / this.coord_rate) + 1 - n_border;
            % Ignore solution at the border of the polynomial
            pid_ceil(pid_ceil < 1) = 1;
            pid_ceil(pid_ceil > size(this.getPolyCoeff, 4)) = size(this.getPolyCoeff, 4);
            
            c_times = this.getCoordTime();
            c_times = c_times - this.time_ref_coord;
            
            poly = permute(poly(:,:,sat_idx, :),[1 3 2 4]);
            
            W_poly = zeros(gps_time.length, 1);
            w = zeros(gps_time.length, 1);
            X_sat = zeros(gps_time.length, n_sat, 3);
            V_sat = zeros(gps_time.length, n_sat, 3);
            n_epoch_old = 0;
            for id = unique([pid_floor; pid_ceil])'
                % find the epochs with the same poly
                p_ids = (pid_floor == id | pid_ceil == id);
                n_epoch = sum(p_ids);
                if n_epoch > 0
                    if n_epoch ~= n_epoch_old
                        t_fct = ones(n_epoch, poly_order + 1);
                        o_mat = repmat(1 : poly_order, n_epoch, 1);
                    end
                    t2 = (t_diff(p_ids) -  c_times(id + n_border))/this.coord_rate;
                    % this is not pre-cached
                    if (n_epoch ~= n_epoch_old) || (sum(abs(t2 - t2_old)) ~= 0)
                        t_fct(:,2) = t2;
                        for o = 3 : poly_order + 1
                            t_fct(:, o) = t_fct(:, o - 1) .* t2;
                        end
                        w = 1 ./ t_fct(:,2) .^ 2;
                        w(t_fct(:, 2) == 0, 1) = 1;
                        w = repmat(w, 1, size(poly, 2), size(poly, 3));
                        n_epoch_old = n_epoch;
                        t2_old = t2;
                    end
                    W_poly(p_ids, 1) = W_poly(p_ids, 1) + w(:,1,1);
                    
                    X_sat(p_ids, :,:) = X_sat(p_ids, :,:) + reshape(t_fct * reshape(poly(:,:,:, id), poly_order + 1, 3 * n_sat), n_epoch, n_sat, 3) .* w;
                    
                    if nargout > 1
                        V_sat(p_ids, :,:) = V_sat(p_ids, :,:) + (reshape((t_fct(:, 1 : poly_order) .* o_mat) * reshape(poly(2 : end, :, :, id), poly_order, 3 * n_sat), n_epoch, n_sat, 3) / this.coord_rate) .* w;
                    end
                end
            end
            %W_poly(W_poly < 1) = nan; % If the polynomial is not stable, do not compute the orbit
            X_sat = X_sat ./ repmat(W_poly, 1, n_sat, 3);
            V_sat = V_sat ./ repmat(W_poly, 1, n_sat, 3);
            
            if size(X_sat,2)==1
                X_sat = squeeze(X_sat);
                V_sat = squeeze(V_sat);
                if size(X_sat,2) == 1
                    X_sat = X_sat';
                    V_sat = V_sat';
                end
            end
        end
        
        function [X_sat, V_sat, A_pol_eval, c_idx] = coordInterpolate11(this, t, sat)
            % SYNTAX:
            %   [X_sat] = Eph_Tab.polInterpolate11(t, sat)
            %
            % INPUT:
            %    t = vector of times where to interpolate
            %    sat = satellite to be interpolated (optional)
            % OUTPUT:
            %
            % DESCRIPTION: interpolate coordinates of staellites expressed with a Lagrange interpolator of degree 11
            n_sat = size(this.coord, 2);
            if nargin <3
                sat_idx = ones(n_sat, 1) > 0;
            else
                sat_idx = sat;
            end
            
            n_sat = length(sat_idx);
            nt = t.length();
            %c_idx=round(t_fd/this.coord_rate)+this.start_time_idx;%coefficient set  index
            poly = this.getPolyCoeff;
            c_idx = round((t - this.time_ref_coord) / this.coord_rate) + 1 - ((size(this.coord, 1) - size(poly, 4)) / 2);
            
            c_idx(c_idx < 1) = 1;
            c_idx(c_idx > size(poly,4)) = size(poly, 4);
            
            c_times = this.getCoordTime();
            % convert to difference from 1st time of the tabulated ephemerids (precise enough in the span of few days and faster that calaling method inside the loop)
            t_diff = t - this.time_ref_coord;
            c_times = c_times - this.time_ref_coord;
            %l_idx=idx-5;
            %u_id=idx+10;
            
            X_sat = zeros(nt,n_sat,3);
            A_pol_eval = zeros(nt,11);
            V_sat = zeros(nt,n_sat,3);
            un_idx = unique(c_idx)';
            for id = un_idx
                t_idx = c_idx == id;
                times = t_diff(t_idx);
                t_fct =  (times -  c_times(id + ((size(this.coord, 1) - size(poly, 4)) / 2)))/this.coord_rate;
                
                %%%% compute position
                t_fct2 = t_fct .* t_fct;
                t_fct3 = t_fct2 .* t_fct;
                t_fct4 = t_fct3 .* t_fct;
                t_fct5 = t_fct4 .* t_fct;
                t_fct6 = t_fct5 .* t_fct;
                t_fct7 = t_fct6 .* t_fct;
                t_fct8 = t_fct7 .* t_fct;
                t_fct9 = t_fct8 .* t_fct;
                t_fct10 = t_fct9 .* t_fct;
                eval_vec = [ones(size(t_fct)) ...
                    t_fct ...
                    t_fct2 ...
                    t_fct3 ...
                    t_fct4 ...
                    t_fct5 ...
                    t_fct6 ...
                    t_fct7 ...
                    t_fct8 ...
                    t_fct9 ...
                    t_fct10];
                X_sat(t_idx,:,:) = reshape(eval_vec*reshape(permute(poly(:,:,sat_idx,id),[1 3 2 4]),11,3*n_sat),sum(t_idx),n_sat,3);
                A_pol_eval(t_idx,:) = A_pol_eval(t_idx,:) + eval_vec;
                %%% compute velocity
                eval_vec = [ ...
                    ones(size(t_fct))  ...
                    2*t_fct  ...
                    3*t_fct2 ...
                    4*t_fct3 ...
                    5*t_fct4 ...
                    6*t_fct5 ...
                    7*t_fct6 ...
                    8*t_fct7 ...
                    9*t_fct8 ...
                    10*t_fct9];
                V_sat(t_idx,:,:) = reshape(eval_vec*reshape(permute(poly(2:end,:,sat_idx,id),[1 3 2 4]),10,3*n_sat),sum(t_idx),n_sat,3)/this.coord_rate;
            end
            if size(X_sat,2)==1
                X_sat = squeeze(X_sat);
                V_sat = squeeze(V_sat);
                if size(X_sat,2) ==1
                    X_sat = X_sat';
                    V_sat = V_sat';
                end
            end
        end
        
        function [sun_ECEF, moon_ECEF ] = sunMoonInterpolate(this, t, no_moon)
            % SYNTAX:
            %   [X_sat]=Eph_Tab.sunInterpolate(t,sat)
            %
            % INPUT:
            %    time = vector of times where to interpolate
            %    no_mmon = do not compute moon postion (default false)
            % OUTPUT:
            %
            % DESCRIPTION: interpolate sun and moon positions
            if isempty(this.X_moon) || isempty(this.X_sun)
                this.tabulateSunMoonPos();
            end
            
            if isempty(this.sun_pol_coeff)
                this.computeSMPolyCoeff();
            end
            
            if nargin < 3
                moon = true;
            else
                moon = not(no_moon);
            end
            %c_idx=round(t_fd/this.coord_rate)+this.start_time_idx;%coefficient set  index
            
            c_idx = round((t - this.time_ref_coord) / this.coord_rate) - 4;
            
            c_idx(c_idx<1) = 1;
            c_idx(c_idx > size(this.X_sun,1)-10) = size(this.X_sun,1)-10;
            
            c_times = this.getCoordTime();
            
            
            %l_idx=idx-5;
            %u_id=idx+10;
            nt = t.length();
            sun_ECEF=zeros(nt,3);
            if moon
                moon_ECEF=zeros(nt,3);
            end
            
            % convert to difference from 1st time of the tabulated ephemerids (precise enough in the span of few days and faster that calaling method inside the loop)
            t = t - this.time_ref_coord;
            c_times = c_times - this.time_ref_coord;
            
            un_idx=unique(c_idx)';
            for idx=un_idx
                t_idx=c_idx==idx;
                times=t(t_idx);
                %t_fct=((times-this.time(5+idx)))';%time from coefficient time
                t_fct =  (times -  c_times(idx+5))/this.coord_rate; %
                %%%% compute position
                t_fct2 = t_fct .* t_fct;
                t_fct3 = t_fct2 .* t_fct;
                t_fct4 = t_fct3 .* t_fct;
                t_fct5 = t_fct4 .* t_fct;
                t_fct6 = t_fct5 .* t_fct;
                t_fct7 = t_fct6 .* t_fct;
                t_fct8 = t_fct7 .* t_fct;
                t_fct9 = t_fct8 .* t_fct;
                t_fct10 = t_fct9 .* t_fct;
                eval_vec = [ones(size(t_fct)) ...
                    t_fct ...
                    t_fct2 ...
                    t_fct3 ...
                    t_fct4 ...
                    t_fct5 ...
                    t_fct6 ...
                    t_fct7 ...
                    t_fct8 ...
                    t_fct9 ...
                    t_fct10];
                sun_ECEF(t_idx,:) = eval_vec*reshape(this.sun_pol_coeff(:,:,idx),11,3);
                if moon
                    moon_ECEF(t_idx,:) = eval_vec*reshape(this.moon_pol_coeff(:,:,idx),11,3);
                end
            end
        end
        
        function [sun_ECEF, moon_ECEF] = computeSunMoonPos(this, time, no_moon)
            % SYNTAX:
            %   this.computeSunMoonPos(p_time)
            %
            % INPUT:
            %    time = Gps_Time [n_epoch x 1]
            %    no_moon = do not compute moon (Boolena deafult false)
            % OUTPUT:
            % sun_ECEF  : sun  coordinate Earth Centered Earth Fixed [n_epoch x 3]
            % moon_ECEF : moon coordinate Earth Centered Earth Fixed [n_epoch x 3]
            % DESCRIPTION: Compute sun and moon psitions at the time
            % desidered time
            
            global iephem km ephname inutate psicor epscor ob2000
            %time = GPS_Time((p_time(1))/86400+GPS_Time.GPS_ZERO);
            if nargin < 3
                moon = true;
            else
                moon = not(no_moon);
            end
            
            sun_id = 11; moon_id = 10; earth_id = 3;
            
            readleap; iephem = 1; ephname = 'de436.bin'; km = 1; inutate = 1; ob2000 = 0.0d0;
            
            tmatrix = j2000_icrs(1);
            
            setmod(2);
            % setdt(3020092e-7);
            setdt(5.877122033683494);
            xp = 171209e-6; yp = 414328e-6;
            
            go_dir = Core.getLocalStorageDir();
            
            %if the binary JPL ephemeris file is not available, try to copy from reference, or generate it
            if (exist(fullfile(go_dir, ephname),'file') ~= 2)
                success = 0;
                if (exist(fullfile('..', 'data', 'reference', 'JPL', ephname),'file') == 2)
                    success = copyfile(fullfile('..', 'data', 'reference', 'JPL', ephname), fullfile(go_dir, ephname), 'f');
                end
                if ~success
                    fprintf('Warning: file "de436.bin" not found in at %s\n         ... generating a new "de436.bin" file\n',fullfile(go_dir, 'de436.bin'));
                    fprintf('         (this procedure may take a while, but it will be done only once on each installation):\n')
                    fprintf('-------------------------------------------------------------------\n\n')
                    asc2eph(436, {'ascp01950.436', 'ascp02050.436'}, fullfile(go_dir, 'de436.bin'));
                    fprintf('-------------------------------------------------------------------\n\n')
                end
            end
            
            sun_ECEF = zeros(time.length(), 3);
            moon_ECEF = zeros(time.length(), 3);
            time = time.getCopy;
            time.toUtc;
            
            jd_utc = time.getJD;
            jd_tdb = time.getJDTDB; % UTC to TDB
            
            % precise celestial pole (disabled)
            %[psicor, epscor] = celpol(jd_tdb, 1, 0.0d0, 0.0d0);
            psicor = 0;
            epscor = 0;
            
            for e = 1 : time.length()
                % compute the Sun position (ICRS coordinates)
                rrd = jplephem(jd_tdb(e), sun_id, earth_id);
                sun_ECI = rrd(1:3);
                sun_ECI = tmatrix * sun_ECI;
                
                % Sun ICRS coordinates to ITRS coordinates
                deltat = getdt;
                jdut1 = jd_utc(e) - deltat;
                tjdh = floor(jdut1); tjdl = jdut1 - tjdh;
                sun_ECEF(e,:) = celter(tjdh, tjdl, xp, yp, sun_ECI)*1e3;
                
                if moon
                    % compute the Moon position (ICRS coordinates)
                    rrd = jplephem(jd_tdb(e), moon_id, earth_id);
                    moon_ECI = rrd(1:3);
                    moon_ECI = tmatrix * moon_ECI;
                    
                    % Moon ICRS coordinates to ITRS coordinates
                    deltat = getdt;
                    jdut1 = jd_utc(e) - deltat;
                    tjdh = floor(jdut1); tjdl = jdut1 - tjdh;
                    moon_ECEF(e,:) = celter(tjdh, tjdl, xp, yp, moon_ECI)*1e3;
                end
            end
        end
        
        function tabulateSunMoonPos(this)
            % SYNTAX:
            %   this.computeSunMoonPos(p_time)
            %
            % INPUT:
            %    p_time = Gps_Time [n_epoch x 1]
            % OUTPUT:
            % DESCRIPTION: Compute sun and moon positions at coordinates time and
            % store them in the object (Overwrite previous data)
            
            %his.t_sun = p_time;
            [this.X_sun , this.X_moon] = this.computeSunMoonPos(this.getCoordTime());
            this.computeSMPolyCoeff();
        end
                
        function loadAntPCO(this)
            % Loading antenna's phase center offsets
            % for code solution
            %
            % SYNTAX
            %   this.loadAntPCO
            
            Core.getLogger.addMessage(Core.getLogger.indent(sprintf('Reading antenna PCV/PCO of %s constellation/s', this.cc.getActiveSysChar)));
            atx = Core.getAntennaManager();
            ant = atx.getAntenna('', this.cc.getAntennaId(), this.time_ref_coord, 100);
            
            this.ant_pco1 = zeros(1, this.cc.getNumSat(), 3);
            if isempty(this.avail)
                this.avail = zeros(this.cc.getNumSat(), 1);
            end
            for sat = 1 : size(ant(:),1)
                if ~ant(sat).isEmpty()
                    if Core.getConstellationCollector.getSysPrn(sat) ~= 'I'
                        this.ant_pco1(:,sat,:) = ant(sat).getPCO([ant(sat).f_code(1) '01'])' * 1e-3;
                    else
                        this.ant_pco1(:,sat,:) = ant(sat).getPCO([ant(sat).f_code(1) '05'])' * 1e-3;
                    end
                else
                    this.avail(sat) = 0;
                end
            end
        end
        
        function writeSP3(this, f_name, prec)
            % SYNTAX:
            %   eph_tab.writeSP3(f_name, prec)
            %
            % INPUT:
            %   f_name       = file name of the sp3 file to be written
            %   prec        = precision (cm) of satellite orbit for all
            %   satellites (default 100)
            %
            %
            % DESCRIPTION:
            %   Write the current satellite postions and clocks bias into a sp3
            %   file
            if nargin<3
                prec=99;
            end
            %%% check if clock rate and coord rate are compatible
            rate_ratio=this.coord_rate/this.clock_rate;
            if abs(rate_ratio-round(rate_ratio)) > 0.00000001
                Core.getLogger.addWarning(sprintf('Incompatible coord rate (%s) and clock rate (%s) , sp3 not produced',this.coord_rate,this.clock_rate))
                return
            end
            %%% check if sun and moon positions have been computed
            if isempty(this.X_sun) || this.X_sun(1,1)==0
                this.sun_moon_pos();
            end
            %%% compute center of mass position (X_sat - PCO)
            switch_back = false;
            if this.poly_type == 1
                this.toCOM();
                switch_back = true;
            end
            %%% write to file
            rate_ratio = round(rate_ratio);
            fid = fopen(f_name,'Wb');
            this.writeHeader(fid, prec);
            
            for i=1:length(this.coord)
                this.writeEpoch(fid,[squeeze(this.coord(i,:,:)/1000) this.clock((i-1)/rate_ratio+1,:)'*1000000],i); %% convert coord in km and clock in microsecodns
            end
            fprintf(fid,'EOF\n');
            fclose(fid);
            if switch_back
                this.toAPC();
            end
        end
        
        function writeHeader(this, fid, prec)
            
            if nargin<3
                %%% unknown precision
                prec=99;
            end
            %prec = num2str(prec);
            time=this.time_ref_coord.getCopy();
            str_time = time.toString();
            year = str2num(str_time(1:4));
            month = str2num(str_time(6:7));
            day = str2num(str_time(9:10));
            hour = str2num(str_time(12:13));
            minute = str2num(str_time(15:16));
            second = str2num(str_time(18:27));
            week = time.getGpsWeek();
            sow = time.getGpsTime()-week*7*86400;
            mjd = jd2mjd(cal2jd(year,month,day));
            d_frac = hour/24+minute/24*60+second/86400;
            step = this.coord_rate;
            num_epoch = length(this.time_ref_coord);
            cc = this.cc;
            fprintf(fid,'#cP%4i %2i %2i %2i %2i %11.8f %7i d+D   IGS14 CNV GReD\n',year,month,day,hour,minute,second,num_epoch);
            fprintf(fid,'## %4i %15.8f %14.8f %5i %15.13f\n',week,sow,step,mjd,d_frac);
            
            sats = [];
            pre = [];
            ids = cc.prn;
            for i = 1:length(ids)
                sats=[sats, strrep(sprintf('%s%2i', cc.system(i), ids(i)), ' ', '0')];
                pre=[pre, sprintf('%3i', prec)];
            end
            n_row=ceil(length(sats)/51);
            rows=cell(5,1);
            rows(:)={repmat('  0',1,17)};
            pres=cell(5,1);
            pres(:)={repmat('  0',1,17)};
            for i =1:n_row
                rows{i}=sats((i-1)*51+1:min(length(sats),i*51));
                pres{i}=pre((i-1)*51+1:min(length(pre),i*51));
            end
            last_row_length=length((i-1)*51+1:length(sats));
            rows{n_row}=[rows{n_row} repmat('  0',1,(51-last_row_length)/3)];
            pres{n_row}=[pres{n_row} repmat('  0',1,(51-last_row_length)/3)];
            
            fprintf(fid,'+   %2i   %s\n',sum(cc.n_sat),rows{1});
            for i=2:length(rows)
                fprintf(fid,'+        %s\n',rows{i});
            end
            for i=1:length(rows)
                fprintf(fid,'++       %s\n',pres{i});
            end
            fprintf(fid,'%%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n');
            fprintf(fid,'%%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n');
            fprintf(fid,'%%f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n');
            fprintf(fid,'%%f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n');
            fprintf(fid,'%%i    0    0    0    0      0      0      0      0         0\n');
            fprintf(fid,'%%i    0    0    0    0      0      0      0      0         0\n');
            fprintf(fid,'/* Produced using goGPS                                     \n');
            fprintf(fid,'/*                 Non                                      \n');
            fprintf(fid,'/*                     Optional                             \n');
            fprintf(fid,'/*                              Lines                       \n');
        end
        
        function writeEpoch(this,fid,XYZT,epoch)
            t = this.time_ref_coord.getCopy();
            t.addIntSeconds((epoch) * 900);
            cc = this.cc;
            str_time = t.toString();
            year = str2num(str_time(1:4));
            month = str2num(str_time(6:7));
            day = str2num(str_time(9:10));
            hour = str2num(str_time(12:13));
            minute = str2num(str_time(15:16));
            second = str2num(str_time(18:27));
            fprintf(fid,'*  %4i %2i %2i %2i %2i %11.8f\n',year,month,day,hour,minute,second);
            for i = 1:size(XYZT,1)
                fprintf(fid,'P%s%14.6f%14.6f%14.6f%14.6f\n',strrep(sprintf('%s%2i', cc.system(i), cc.prn(i)), ' ', '0'),XYZT(i,1),XYZT(i,2),XYZT(i,3),XYZT(i,4));
            end
            
        end
        
        function sys_c = getAvailableSys(this)
            % get the available system stored into the object
            % SYNTAX: sys_c = this.getAvailableSys()
            
            % Select only the systems present in the file
            sys_c = this.cc.getAvailableSys();
        end
    end
    
    % ==================================================================================================================================================
    %% STATIC FUNCTIONS used as utilities
    % ==================================================================================================================================================
    methods (Static, Access = public)
        
        function [az, el, sat_coo] = getAzEl(coo, time, sat)
            % Get the azimuth and elevation of a satellite "sat" e.g. G21
            % 
            % INPUT
            %   coo         Coordinates of the ground point
            %   time        GPS_Time object
            %   sat         satellite (eg. 'R02') / go_id (e.g. 34)
            %
            % OUTPUT
            %   az, el      azimuth and elevation [deg]
            %   sat_coo     coordinates of the satellite (Coordinate object)
            %
            % SYNTAX
            %   [az, el, sat_coo] = Core_Sky.getAzEl(coo, gps_time, sat)
            %
            % EXAMPLE 
            %   [az, el] = Core_Sky.getAzEl(Coordinates.fromGeodetic(0,0,0), GPS_Time('2020-10-14 12:14'), 'G22');
            
            core = Core.getCurrentCore;
            old_cc = unique(core.getConstellationCollector.getActiveSysChar);
            sky = core.sky;
            if isempty(sky)
                sky = Core_Sky();
                core.sky = sky;
            end
            cc = sky.cc;
            if isempty(sat)
                sat = cc.index;
            end
            if isnumeric(sat)
                sat_name = cc.getSatName(sat);
            else
                sat_name = sat;
            end
            sys_char_list = unique(sat_name(:,1)');
            cc.setActive(sys_char_list);
            if ~strcmp(sys_char_list, old_cc)
                core.state.cc = cc;
                fw = File_Wizard();                
                if (time.last < (GPS_Time.now.addIntSeconds(-86400*14))); 
                    fw.setCurCenter('code_mgex_aiub', 'final');
                else
                    fw.setCurCenter('code_predicted', 'ultra');
                end
                fw.conjureNavFiles(time.first, time.last);
            else
                core.state.cc = cc;            
                if isempty(core.state.eph_name)
                    fw = File_Wizard();
                    fw.conjureNavFiles(time.first, time.last);
                end
            end
            
            lim = time.first.getCopy;
            lim.append(time.last);
            flag_no_clock = true;
            core.initSkySession(lim, flag_no_clock);
            if nargout == 3
                [az, el, sat_coo] = sky.getAzimuthElevation(coo, time, sat);
            else
                [az, el] = sky.getAzimuthElevation(coo, time, sat);
            end
            
        end
        
        function prn_num = prnName2Num(prn_name)
            % Convert a 4 char name into a numeric value (float)
            % SYNTAX:
            %   marker_num = markerName2Num(marker_name);
            
            prn_num = prn_name(:,1:3) * [2^16 2^8 1]';
        end
        
        function prn_name = prnNum2Name(prn_num)
            % Convert a numeric value (float) of a station into a 4 char marker
            % SYNTAX:
            %   marker_name = markerNum2Name(marker_num)
            prn_name = char(zeros(numel(prn_num), 3));
            prn_name(:,1) = char(floor(prn_num / 2^16));
            prn_num = prn_num - prn_name(:,1) * 2^16;
            prn_name(:,2) = char(floor(prn_num / 2^8));
            prn_num = prn_num - prn_name(:,2) * 2^8;
            prn_name(:,3) = char(prn_num);
        end
        
        function [eph, iono] = loadNavParameters(file_nav, cc)
            % SYNTAX:
            %   [eph, iono] = getNavParameters(file_nav, cc);
            %
            % INPUT:
            %   file_nav = RINEX navigation file
            %   cc = Constellation_Collector object, contains the satus of the satellite systems in use
            %
            % OUTPUT:
            %   Eph = matrix containing 33 navigation parameters for each satellite
            %   iono = matrix containing ionosphere parameters
            %
            % DESCRIPTION:
            %   Parse a RINEX navigation file.
            
            %  Partially based on RINEXE.M (EASY suite) by Kai Borre
            
            % ioparam = 0;
            eph = [];
            iono = zeros(8,1);
            
            
            %%
            % open RINEX observation file
            fid = fopen(file_nav,'rt');
            if fid < 0
                Core.getLogger.addWarning(sprintf('File not found "%s", unable to import navigation parameters', file_nav))
            else
                txt = fread(fid,'*char')';
                % try to see if carriage return is present in the file (Windows stupid standard)
                % On Windows file lines ends with char(13) char(10)
                % instead of just using char(10)
                if ~isempty(find(txt(1:min(1000,numel(txt))) == 13, 1, 'first'))
                    has_cr = true;  % The file has carriage return - I hate you Bill!
                else
                    has_cr = false;  % The file is UNIX standard
                end
                % txt = txt(txt ~= 13);  % remove carriage return - I hate you Bill!
                fclose(fid);
                
                % get new line separators
                nl = regexp(txt, '\n')';
                if nl(end) <  (numel(txt) - double(has_cr))
                    nl = [nl; numel(txt)];
                end
                lim = [[1; nl(1 : end - 1) + 1] (nl - 1 - double(has_cr))];
                lim = [lim lim(:,2) - lim(:,1)];
                while lim(end,3) < 3
                    lim(end,:) = [];
                end
                
                % removing empty lines at end of file
                while (lim(end,1) - lim(end-1,1))  < 2
                    lim(end,:) = [];
                end
                
                eoh = 0;
                flag_eoh = false;
                while eoh < size(lim, 1) && flag_eoh == false
                    eoh = eoh + 1;
                    flag_eoh = strcmp(txt((lim(eoh,1) + 60) : min(lim(eoh,1) + 72, lim(eoh, 2))), 'END OF HEADER');
                end
                
                % Reading Header
                head_field{1} = 'RINEX VERSION / TYPE';                  %  1
                head_field{2} = 'PGM / RUN BY / DATE';                   %  2
                head_field{3} = 'LEAP SECONDS';                          %  3
                head_field{4} = 'ION ALPHA';                             %  4
                head_field{5} = 'ION BETA';                              %  5
                
                % parsing ------------------------------------------------------------------------------------------------------------------------------------------
                
                % retriving the kind of header information is contained on each line
                line2head = zeros(eoh, 1);
                l = 0;
                while l < eoh
                    l = l + 1;
                    %DEBUG: txt((lim(l,1) + 60) : lim(l,2))
                    tmp = find(strcmp(strtrim(txt((lim(l,1) + 60) : lim(l,2))), head_field));
                    if ~isempty(tmp)
                        % if the field have been recognized (it's not a comment)
                        line2head(l) = tmp;
                    end
                end
                
                % read RINEX type 3 or 2 ---------------------------------------------------------------------------------------------------------------------------
                
                l = find(line2head == 1);
                type_found = ~isempty(l);
                
                if type_found
                    dataset = textscan(txt(lim(1,1):lim(1,2)), '%f%c%18c%c');
                else
                    throw(MException('VerifyINPUTInvalidNavigationalFile', 'This navigational RINEX does not contain orbits'));
                end
                this.rin_type = dataset{1};
                this.rinex_ss = dataset{4};
                
                if dataset{2} == 'N'
                    % Is a navigational file
                else
                    throw(MException('VerifyINPUTInvalidNavigationalFile', 'This navigational RINEX does not contain orbits'));
                end
                
                % Read iono parameters (if found):
                
                [~, l_iono] = intersect(line2head, [4,5]);
                iono_found = numel(l_iono) == 2;
                iono_loaded = false;
                if iono_found
                    data = textscan(txt(lim(l_iono(1),1) + (2 : 49)), '%f%f%f%f');
                    if ~isempty(data{4})
                        iono(1) = data{1};
                        iono(2) = data{2};
                        iono(3) = data{3};
                        iono(4) = data{4};
                    end
                    data = textscan(txt(lim(l_iono(2),1) + (2 : 49)), '%f%f%f%f');
                    if ~isempty(data{4})
                        iono(5) = data{1};
                        iono(6) = data{2};
                        iono(7) = data{3};
                        iono(8) = data{4};
                    end
                    iono_loaded = true;
                end
                
                if this.rin_type < 3 % at the moment the new reader support only RINEX 3 broadcast ephemeris
                    [eph, iono] = RINEX_get_nav(file_nav, cc); % Old implementation slower but support RINEX 2
                else
                    eph = [];
                    for sys_c = cc.getActiveSysChar()
                        
                        id_ss = find(txt(lim(eoh:end,1)) == sys_c) + eoh - 1;
                        n_epo = numel(id_ss);
                        
                        switch sys_c
                            case 'G'
                                sys_index = cc.getGPS().getFirstId();
                                nppl = [4 4 4 4 4 4 4 2]; % number of parameters per line
                            case 'R'
                                sys_index = cc.getGLONASS().getFirstId();
                                nppl = [4 4 4 4];         % number of parameters per line
                            case 'E'
                                sys_index = cc.getGalileo().getFirstId();
                                full_line = median(lim(id_ss + 6 ,3)) > 65; % detect if all the lines have 79 chars (empty fields are filled with spaces)
                                nppl = [4 4 4 4 4 (3 + full_line) 4 1]; % number of parameters per line
                            case 'J'
                                sys_index = cc.getQZSS().getFirstId();
                                nppl = [4 4 4 4 4 4 4 2]; % number of parameters per line
                            case 'C'
                                sys_index = cc.getBeiDou().getFirstId();
                                nppl = [4 4 4 4 4 4 4 2]; % number of parameters per line
                            case 'I'
                                sys_index = cc.getIRNSS().getFirstId();
                                full_line = median(lim(id_ss + 6 ,3)) > 65; % detect if all the lines have 79 chars (empty fields are filled with spaces)
                                nppl = [4 4 4 4 4 (3 + full_line) (3 + full_line) 1]; % number of parameters per line
                            case 'S'
                                sys_index = cc.getSBAS().getFirstId();
                                nppl = [4 4 4 4];         % number of parameters per line
                        end
                        par_offset = [4 23 42 61]; % character offset for reading a parameter
                        lin_offset = [0 cumsum((nppl(1:end-1) * 19 + 5 + has_cr))]; % character offset for reading on a certain line
                        
                        % Function to extract a parameter from the broadcast info table
                        getParStr = @(r,c) txt(repmat(lim(id_ss,1),1,19) + par_offset(c) + lin_offset(r) + repmat(0:18, n_epo, 1));
                        getParNum = @(r,c) str2num(txt(repmat(lim(id_ss,1),1,19) + par_offset(c) + lin_offset(r) + repmat(0:18, n_epo, 1)));
                        
                        % Epochs
                        eph_ss = zeros(n_epo, 33);
                        eph_ss(:,  1) = str2num(txt(repmat(lim(id_ss,1), 1, 2) + repmat([1 2], length(id_ss), 1)));
                        
                        date = cell2mat(textscan(getParStr(1,1)','%4f %2f %2f %2f %2f %2f'));
                        time = GPS_Time(date, [], iif(sys_c == 'R', false, true));
                        
                        % Other parameters
                        if ismember(sys_c, 'RS')
                            eph_ss(:,  2) = -getParNum(1,2); % TauN
                            eph_ss(:,  3) = getParNum(1,3); % GammaN
                            eph_ss(:,  4) = getParNum(1,4); % tk
                            
                            eph_ss(:,  5) = 1e3 * getParNum(2,1); % X
                            eph_ss(:,  8) = 1e3 * getParNum(2,2); % Xv
                            eph_ss(:, 11) = 1e3 * getParNum(2,3); % Xa
                            eph_ss(:, 27) = getParNum(2,4); % Bn
                            
                            eph_ss(:,  6) = 1e3 * getParNum(3,1); % Y
                            eph_ss(:,  9) = 1e3 * getParNum(3,2); % Yv
                            eph_ss(:, 12) = 1e3 * getParNum(3,3); % Ya
                            eph_ss(:, 15) = getParNum(3,4); % freq_num
                            
                            eph_ss(:,  7) = 1e3 * getParNum(4,1); % Z
                            eph_ss(:, 10) = 1e3 * getParNum(4,2); % Zv
                            eph_ss(:, 13) = 1e3 * getParNum(4,3); % Za
                            eph_ss(:, 14) = getParNum(4,4); % E
                            
                            [week_toe, toe] = time.getGpsWeek;
                            eph_ss(:, 18) = toe;
                            eph_ss(:, 24) = week_toe;
                            eph_ss(:, 32) = double(week_toe) * 7 * 86400 + toe;
                            
                            eph_ss(:, 30) = eph_ss(:, 1) + (sys_index - 1); % go_id
                            eph_ss(:, 31) = int8(sys_c);
                        else % for GEJCI
                            eph_ss(:, 19) = getParNum(1,2); % af0
                            eph_ss(:, 20) = getParNum(1,3); % af1
                            eph_ss(:,  2) = getParNum(1,4); % af2
                            
                            eph_ss(:, 22) = getParNum(2,1); % IODE
                            eph_ss(:, 11) = getParNum(2,2); % crs
                            eph_ss(:,  5) = getParNum(2,3); % deltan
                            eph_ss(:,  3) = getParNum(2,4); % M0
                            
                            eph_ss(:,  8) = getParNum(3,1); % cuc
                            eph_ss(:,  6) = getParNum(3,2); % ecc
                            eph_ss(:,  9) = getParNum(3,3); % cus
                            eph_ss(:,  4) = getParNum(3,4); % roota
                            
                            eph_ss(:, 18) = getParNum(4,1); % toe
                            eph_ss(:, 14) = getParNum(4,2); % cic
                            eph_ss(:, 16) = getParNum(4,3); % Omega0
                            eph_ss(:, 15) = getParNum(4,4); % cis
                            
                            eph_ss(:, 12) = getParNum(5,1); % i0
                            eph_ss(:, 10) = getParNum(5,2); % crc
                            eph_ss(:,  7) = getParNum(5,3); % omega
                            eph_ss(:, 17) = getParNum(5,4); % Omegadot
                            
                            eph_ss(:, 13) = getParNum(6,1); % idot
                            eph_ss(:, 23) = iif(isempty(getParNum(6,2)),0,getParNum(6,2)); % code_on_L2
                            if (sys_c == 'C') % Beidou week have an offset of 1356 weeks
                                eph_ss(:, 24) = GPS_Time.GPS_BDS_WEEK0 + getParNum(6,3); % weekno
                            else
                                eph_ss(:, 24) = getParNum(6,3); % weekno
                            end
                            if ismember(sys_c, 'GJC') % present only for G,J,C constellations
                                eph_ss(:, 25) = getParNum(6,4); % L2flag
                            end
                            
                            eph_ss(:, 26) = getParNum(7,1); % svaccur
                            eph_ss(:, 27) = getParNum(7,2); % svhealth
                            eph_ss(:, 28) = getParNum(7,3); % tgd
                            
                            %eph_ss(:, xx) = getParNum(8,1); % tom
                            if ismember(sys_c, 'GJC') % present only for G,J,C constellations
                                valid_fit_int = any(getParStr(8, 2)' - 32);
                                eph_ss(valid_fit_int, 29) = getParNum(8, 2); % fit_int
                            end
                            
                            % Other parameter to stor in eph
                            eph_ss(:, 30) = eph_ss(:, 1) + (sys_index - 1); % go_id
                            eph_ss(:, 31) = int8(sys_c);
                            
                            [week, toc] = time.getGpsWeek;
                            eph_ss(:, 21) = toc;
                            eph_ss(:, 32) = double(week)*7*86400 + eph_ss(:, 18);
                            eph_ss(:, 33) = time.getGpsTime();
                            
                            if ismember(sys_c, 'G') % present only for G constellation
                                iodc = getParNum(7,4); % IODC
                                
                                time_thr = 0;
                                iod_check = (abs(eph_ss(:, 22) - iodc) > time_thr);
                                sat_ko = unique(eph_ss(iod_check, 1));
                                log = Core.getLogger;
                                %cm = log.getColorMode();
                                %log.setColorMode(0);
                                if not(isempty(sat_ko))
                                    log.addWarning(sprintf('IODE - IODC of sat %s are different!\nPossible problematic broadcast orbits found for "%s"\nignoring those satellites', sprintf('G%02d ', sat_ko), File_Name_Processor.getFileName(file_nav)));
                                    eph_ss(iod_check, :) = []; % delete non valid ephemeris
                                end
                                %log.setColorMode(cm);
                            end
                        end
                        % Append SS ephemeris
                        eph = [eph eph_ss'];
                    end
                end
            end
        end
        
        % ---------------------------------------------------------------------------
        % Old goGPS functions , integrated with minor modifications as static methods
        %----------------------------------------------------------------------------
                
        function [Eph, iono, flag_return] = loadRinexNav(filename, cc, flag_SP3, iono_model, time, wait_dlg)
            
            % SYNTAX:
            %   [Eph, iono, flag_return] = loadRinexNav(filename, constellations, flag_SP3, iono_model, time, wait_dlg);
            %
            % INPUT:
            %   filename = RINEX navigation file
            %   cc = Constellation_Collector object, contains the satus of the satellite systems in use
            %   flag_SP3 = boolean flag to indicate SP3 availability
            %   wait_dlg = optional handler to waitbar figure (optional)
            %
            % OUTPUT:
            %   Eph = matrix containing 33 navigation parameters for each satellite
            %   iono = vector containing ionosphere parameters
            %   flag_return = notify the parent function that it should return
            %                 (downloaded navigation file still compressed).
            %
            % DESCRIPTION:
            %   Parses RINEX navigation files.
            
            % Check the input arguments
            if (nargin < 6)
                wait_dlg_PresenceFlag = false;
            else
                wait_dlg_PresenceFlag = true;
            end
            
            if (iscell(filename))
                filename = filename{1};
            end
            
            flag_return = 0;
            log = Logger.getInstance();
            state = Core.getCurrentSettings();
            
            %number of satellite slots for enabled constellations
            nSatTot = cc.getNumSat();
            
            %read navigation files
            if (~flag_SP3)
                parse_file(0);
            else
                Eph = zeros(33,nSatTot);
                iono = zeros(8,1);
            end
            
            % Broadcast corrections in DD are currently causing problems (offset in UP) => not using them
            %if Klobuchar ionospheric delay correction is requested but parameters are not available in the navigation file, try to download them
            if ((iono_model == 2 && ~any(iono)) || (flag_SP3 && cc.getGLONASS().isActive()))
                [year, DOY] = time.getDOY();
                
                filename_brdm = ['brdm' num2str(DOY,'%03d') '0.' num2str(mod(year, 100),'%02d') 'p'];
                filename_brdc = ['brdc' num2str(DOY,'%03d') '0.' num2str(mod(year, 100),'%02d') 'n'];
                filename_CGIM = ['CGIM' num2str(DOY,'%03d') '0.' num2str(mod(year, 100),'%02d') 'N'];
                
                pos = find(filename == '/'); if(isempty(pos)), pos = find(filename == '\'); end
                nav_path = filename(1:pos(end));
                
                flag_GLO = flag_SP3 && cc.getGLONASS().isActive();
                
                file_avail = 0;
                if (exist([nav_path filename_brdm],'file') && flag_GLO)
                    filename = [nav_path filename_brdm];
                    file_avail = 1;
                elseif (exist([nav_path filename_CGIM],'file') && ~flag_GLO)
                    filename = [nav_path filename_CGIM];
                    file_avail = 1;
                elseif (exist([nav_path filename_brdc],'file') && ~flag_GLO)
                    filename = [nav_path filename_brdc];
                    file_avail = 1;
                else
                    if (flag_GLO)
                        filename = filename_brdm;
                    else
                        filename = filename_brdc;
                    end
                    [download_successful, compressed] = download_nav(filename, nav_path);
                    filename = [nav_path filename];
                    if (download_successful)
                        file_avail = 1;
                    end
                    if (compressed)
                        flag_return = 1;
                    end
                end
                
                if (file_avail)
                    if (flag_GLO)
                        only_iono = 0;
                    else
                        only_iono = 1;
                    end
                    parse_file(only_iono);
                end
            end
            
            function parse_file(only_iono)
                
                if (wait_dlg_PresenceFlag)
                    waitbar(0.5,wait_dlg,'Reading navigation files...')
                end
                
                Eph_G = []; iono_G = zeros(8,1);
                Eph_R = []; iono_R = zeros(8,1);
                Eph_E = []; iono_E = zeros(8,1);
                Eph_C = []; iono_C = zeros(8,1);
                Eph_J = []; iono_J = zeros(8,1);
                Eph_I = []; iono_I = zeros(8,1);
                
                if (strcmpi(filename(end),'p'))
                    flag_mixed = 1;
                else
                    flag_mixed = 0;
                end
                
                if (cc.getGPS().isActive() || flag_mixed || only_iono)
                    if (exist(filename,'file'))
                        %parse RINEX navigation file (GPS) NOTE: filename expected to
                        %end with 'n' or 'N' (GPS) or with 'p' or 'P' (mixed GNSS)
                        if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename ': ... '])); end
                        % [Eph_G, iono_G] = RINEX_get_nav(filename, cc); % Old implementation slower but support RINEX 2
                        [Eph_G, iono_G] = Core_Sky.loadNavParameters(filename, cc);
                        for sys_c = iif(flag_mixed, cc.SYS_C(cc.active_list), 'G')
                            % Detect and remove satellites with "high" PRN,
                            % usually connected with satellites under testing
                            if ~isempty(Eph_G)
                                id_testing = Eph_G(1,:) > cc.getSys(sys_c).N_SAT & Eph_G(31,:) == sys_c;
                                Eph_G(:, id_testing) = [];
                            else
                                % Core.getLogger.addWarning(sprintf('Ionospheric parameters not found in "%s"', filename));
                            end
                        end
                        if(~only_iono), log.addStatusOk(); end
                    else
                        log.addWarning('GPS navigation file not found. GPS positioning may not work. \n');
                        % cc.deactivateGPS();
                    end
                end
                
                if (cc.getGLONASS().isActive() && ~only_iono)
                    if (exist([filename(1:end-1) 'g'],'file'))
                        %parse RINEX navigation file (GLONASS)
                        if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename(1:end-1) 'g: ... '])); end
                        [Eph_R, iono_R] = Core_Sky.loadNavParameters([filename(1:end-1) 'g'], cc);
                        % Detect and remove satellites with "high" PRN,
                        % usually connected with satellites under testing
                        id_testing = Eph_R(1,:) > GLONASS_SS.N_SAT & Eph_R(31,:) == 'R';
                        Eph_R(:, id_testing) = [];
                        if(~only_iono), log.addStatusOk(); end
                    elseif (~flag_mixed)
                        log.addWarning('GLONASS navigation file not found. GLONASS positioning may not work. \n');
                        % cc.deactivateGLONASS();
                    end
                end
                
                if (cc.getGalileo().isActive() && ~only_iono)
                    if (exist([filename(1:end-1) 'l'],'file'))
                        %parse RINEX navigation file (Galileo)
                        if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename(1:end-1) 'l: ... '])); end
                        [Eph_E, iono_E] = Core_Sky.loadNavParameters([filename(1:end-1) 'l'], cc);
                        % Detect and remove satellites with "high" PRN,
                        % usually connected with satellites under testing
                        id_testing = Eph_E(1,:) > GLONASS_SS.N_SAT & Eph_E(31,:) == 'E';
                        Eph_E(:, id_testing) = [];
                        if(~only_iono), log.addStatusOk(); end
                    elseif (~flag_mixed)
                        log.addWarning('Galileo navigation file not found. Galileo positioning may not work. \n');
                        % cc.deactivateGalileo();
                    end
                end
                
                if (cc.getBeiDou().isActive() && ~only_iono)
                    if (exist([filename(1:end-1) 'c'],'file'))
                        %parse RINEX navigation file (BeiDou)
                        if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename(1:end-1) 'c: ... '])); end
                        [Eph_C, iono_C] = Core_Sky.loadNavParameters([filename(1:end-1) 'c'], cc);
                        % Detect and remove satellites with "high" PRN,
                        % usually connected with satellites under testing
                        id_testing = Eph_C(1,:) > GLONASS_SS.N_SAT & Eph_C(31,:) == 'C';
                        Eph_C(:, id_testing) = [];
                        if(~only_iono), log.addStatusOk(); end
                    elseif (~flag_mixed)
                        log.addWarning('BeiDou navigation file not found. BeiDou positioning may not work. \n');
                        % cc.deactivateBeiDou();
                    end
                end
                
                if (cc.getQZSS().isActive() && ~only_iono)
                    if (exist([filename(1:end-1) 'q'],'file'))
                        %parse RINEX navigation file (QZSS)
                        if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename(1:end-1) 'q: ... '])); end
                        [Eph_J, iono_J] = Core_Sky.loadNavParameters([filename(1:end-1) 'q'], cc);
                        % Detect and remove satellites with "high" PRN,
                        % usually connected with satellites under testing
                        id_testing = Eph_J(1,:) > GLONASS_SS.N_SAT & Eph_J(31,:) == 'J';
                        Eph_J(:, id_testing) = [];
                        if(~only_iono), log.addStatusOk(); end
                    elseif (~flag_mixed)
                        log.addWarning('QZSS navigation file not found. QZSS positioning may not work. \n');
                        % cc.deactivateQZSS();
                    end
                end
                
                if (cc.getIRNSS().isActive() && ~only_iono)
                    if (exist([filename(1:end-1) 'i'],'file'))
                        %parse RINEX navigation file (IRNSS)
                        if(~only_iono), log.addMessage(sprintf('%s',['Reading RINEX file ' filename(1:end-1) 'q: ... '])); end
                        [Eph_I, iono_I] = Core_Sky.loadNavParameters([filename(1:end-1) 'i'], cc);
                        % Detect and remove satellites with "high" PRN,
                        % usually connected with satellites under testing
                        id_testing = Eph_I(1,:) > GLONASS_SS.N_SAT & Eph_I(31,:) == 'I';
                        Eph_I(:, id_testing) = [];
                        if(~only_iono), log.addStatusOk(); end
                    elseif (~flag_mixed)
                        log.addWarning('IRNSS navigation file not found. QZSS positioning may not work. \n');
                        % cc.deactivateIRNSS();
                    end
                end
                
                if (~only_iono)
                    Eph = [Eph_G Eph_R Eph_E Eph_C Eph_J Eph_I];
                end
                
                if (any(iono_G))
                    iono = iono_G;
                elseif (any(iono_R))
                    iono = iono_R;
                elseif (any(iono_E))
                    iono = iono_E;
                elseif (any(iono_C))
                    iono = iono_C;
                elseif (any(iono_J))
                    iono = iono_J;
                elseif (any(iono_I))
                    iono = iono_I;
                else
                    iono = zeros(8,1);
                    if isempty(regexp(filename, '(?<=brdm).*', 'once')) % brdm are broadcast mgex with no iono parameters, iono will be imported from other files
                        log.addWarning(sprintf('Klobuchar ionosphere parameters not found in navigation file\n("%s")\n', filename));
                    end
                end
                
                if (wait_dlg_PresenceFlag)
                    waitbar(1,wait_dlg)
                end
            end
        end
        
        function [satp, satv] = satelliteOrbits(t, Eph, sat, sbas)
            
            % SYNTAX:
            %   [satp, satv] = satelliteOrbits(t, Eph, sat, sbas);
            %
            % INPUT:
            %   t = clock-corrected GPS time
            %   Eph  = ephemeris matrix
            %   sat  = satellite index
            %   sbas = SBAS corrections
            %
            % OUTPUT:
            %   satp = satellite position (X,Y,Z)
            %   satv = satellite velocity
            %
            % DESCRIPTION:
            %   Computation of the satellite position (X,Y,Z) and velocity by means
            %   of its ephemerides.
            
            % the following two line offer an elegant but slow implementation
            %cc = Constellation_Collector('GRECJI');
            %sys_str = cc.getSys(char(Eph(31)));            
            % the following switch is equivalent to the previous two lines but musch faster
            switch char(Eph(31))
                case 'G'
                    sys_str = GPS_SS();
                case 'R'
                    sys_str = GLONASS_SS();
                case 'E'
                    sys_str = Galileo_SS();
                case 'C'
                    sys_str = BeiDou_SS();
                case 'J'
                    sys_str = QZSS_SS();
                case 'I'
                    sys_str = IRNSS_SS();
                case 'S'
                    sys_str = SBAS_SS();
            end
            
            orbital_p = sys_str.ORBITAL_P;
            Omegae_dot = orbital_p.OMEGAE_DOT;
            
            
            
            %consider BeiDou time (BDT) for BeiDou satellites
            if (strcmp(char(Eph(31)),'C'))
                t = t - 14;
            end
            
            %GPS/Galileo/BeiDou/QZSS satellite coordinates computation
            if (~strcmp(char(Eph(31)),'R'))
                
                %get ephemerides
                roota     = Eph(4);
                ecc       = Eph(6);
                omega     = Eph(7);
                cuc       = Eph(8);
                cus       = Eph(9);
                crc       = Eph(10);
                crs       = Eph(11);
                i0        = Eph(12);
                IDOT      = Eph(13);
                cic       = Eph(14);
                cis       = Eph(15);
                Omega0    = Eph(16);
                Omega_dot = Eph(17);
                toe       = Eph(18);
                time_eph  = Eph(32);
                
                %SBAS satellite coordinate corrections
                if (~isempty(sbas))
                    dx_sbas = sbas.dx(sat);
                    dy_sbas = sbas.dy(sat);
                    dz_sbas = sbas.dz(sat);
                else
                    dx_sbas = 0;
                    dy_sbas = 0;
                    dz_sbas = 0;
                end
                
                %-------------------------------------------------------------------------------
                % ALGORITHM FOR THE COMPUTATION OF THE SATELLITE COORDINATES (IS-GPS-200E)
                %-------------------------------------------------------------------------------
                
                %eccentric anomaly
                [Ek, n] = ecc_anomaly(t, Eph);
                
                cr = 6.283185307179600;
                
                A = roota*roota;             %semi-major axis
                tk = check_t(t - time_eph);  %time from the ephemeris reference epoch
                
                fk = atan2(sqrt(1-ecc^2)*sin(Ek), cos(Ek) - ecc);    %true anomaly
                phik = fk + omega;                           %argument of latitude
                phik = rem(phik,cr);
                
                uk = phik                + cuc*cos(2*phik) + cus*sin(2*phik); %corrected argument of latitude
                rk = A*(1 - ecc*cos(Ek)) + crc*cos(2*phik) + crs*sin(2*phik); %corrected radial distance
                ik = i0 + IDOT*tk        + cic*cos(2*phik) + cis*sin(2*phik); %corrected inclination of the orbital plane
                
                %satellite positions in the orbital plane
                x1k = cos(uk)*rk;
                y1k = sin(uk)*rk;
                
                %if GPS/Galileo/QZSS or MEO/IGSO BeiDou satellite
                if (~strcmp(char(Eph(31)),'C') || (strcmp(char(Eph(31)),'C') && Eph(1) > 5))
                    
                    %corrected longitude of the ascending node
                    Omegak = Omega0 + (Omega_dot - Omegae_dot)*tk - Omegae_dot*toe;
                    Omegak = rem(Omegak + cr, cr);
                    
                    %satellite Earth-fixed coordinates (X,Y,Z)
                    xk = x1k*cos(Omegak) - y1k*cos(ik)*sin(Omegak);
                    yk = x1k*sin(Omegak) + y1k*cos(ik)*cos(Omegak);
                    zk = y1k*sin(ik);
                    
                    %apply SBAS corrections (if available)
                    satp = zeros(3,1);
                    satp(1,1) = xk + dx_sbas;
                    satp(2,1) = yk + dy_sbas;
                    satp(3,1) = zk + dz_sbas;
                    
                else %if GEO BeiDou satellite (ranging code number <= 5)
                    
                    %corrected longitude of the ascending node
                    Omegak = Omega0 + Omega_dot*tk - Omegae_dot*toe;
                    Omegak = rem(Omegak + cr, cr);
                    
                    %satellite coordinates (X,Y,Z) in inertial system
                    xgk = x1k*cos(Omegak) - y1k*cos(ik)*sin(Omegak);
                    ygk = x1k*sin(Omegak) + y1k*cos(ik)*cos(Omegak);
                    zgk = y1k*sin(ik);
                    
                    %store inertial coordinates in a vector
                    Xgk = [xgk; ygk; zgk];
                    
                    %rotation matrices from inertial system to CGCS2000
                    Rx = [1        0          0;
                        0 +cosd(-5) +sind(-5);
                        0 -sind(-5) +cosd(-5)];
                    
                    oedt = Omegae_dot*tk;
                    
                    Rz = [+cos(oedt) +sin(oedt) 0;
                        -sin(oedt) +cos(oedt) 0;
                        0           0         1];
                    
                    %apply the rotations
                    Xk = Rz*Rx*Xgk;
                    
                    xk = Xk(1);
                    yk = Xk(2);
                    zk = Xk(3);
                    
                    %store CGCS2000 coordinates
                    satp = zeros(3,1);
                    satp(1,1) = xk;
                    satp(2,1) = yk;
                    satp(3,1) = zk;
                end
                
                %-------------------------------------------------------------------------------
                % ALGORITHM FOR THE COMPUTATION OF THE SATELLITE VELOCITY (as in Remondi,
                % GPS Solutions (2004) 8:181-183 )
                %-------------------------------------------------------------------------------
                if (nargout > 1)
                    Mk_dot = n;
                    Ek_dot = Mk_dot/(1-ecc*cos(Ek));
                    fk_dot = sin(Ek)*Ek_dot*(1+ecc*cos(fk)) / ((1-cos(Ek)*ecc)*sin(fk));
                    phik_dot = fk_dot;
                    uk_dot = phik_dot + 2*(cus*cos(2*phik)-cuc*sin(2*phik))*phik_dot;
                    rk_dot = A*ecc*sin(Ek)*Ek_dot + 2*(crs*cos(2*phik)-crc*sin(2*phik))*phik_dot;
                    ik_dot = IDOT + 2*(cis*cos(2*phik)-cic*sin(2*phik))*phik_dot;
                    Omegak_dot = Omega_dot - Omegae_dot;
                    x1k_dot = rk_dot*cos(uk) - y1k*uk_dot;
                    y1k_dot = rk_dot*sin(uk) + x1k*uk_dot;
                    xk_dot = x1k_dot*cos(Omegak) - y1k_dot*cos(ik)*sin(Omegak) + y1k*sin(ik)*sin(Omegak)*ik_dot - yk*Omegak_dot;
                    yk_dot = x1k_dot*sin(Omegak) + y1k_dot*cos(ik)*cos(Omegak) - y1k*sin(ik)*ik_dot*cos(Omegak) + xk*Omegak_dot;
                    zk_dot = y1k_dot*sin(ik) + y1k*cos(ik)*ik_dot;
                    
                    satv = zeros(3,1);
                    satv(1,1) = xk_dot;
                    satv(2,1) = yk_dot;
                    satv(3,1) = zk_dot;
                end
                
            else %GLONASS satellite coordinates computation (GLONASS-ICD 5.1)
                
                time_eph = Eph(32); %ephemeris reference time
                
                X   = Eph(5);  %satellite X coordinate at ephemeris reference time
                Y   = Eph(6);  %satellite Y coordinate at ephemeris reference time
                Z   = Eph(7);  %satellite Z coordinate at ephemeris reference time
                
                Xv  = Eph(8);  %satellite velocity along X at ephemeris reference time
                Yv  = Eph(9);  %satellite velocity along Y at ephemeris reference time
                Zv  = Eph(10); %satellite velocity along Z at ephemeris reference time
                
                Xa  = Eph(11); %acceleration due to lunar-solar gravitational perturbation along X at ephemeris reference time
                Ya  = Eph(12); %acceleration due to lunar-solar gravitational perturbation along Y at ephemeris reference time
                Za  = Eph(13); %acceleration due to lunar-solar gravitational perturbation along Z at ephemeris reference time
                %NOTE:  Xa,Ya,Za are considered constant within the integration interval (i.e. toe ?}15 minutes)
                
                %integration step
                int_step = 60; %[s]
                
                %time from the ephemeris reference epoch
                tk = check_t(t - time_eph);
                
                %number of iterations on "full" steps
                n = floor(abs(tk/int_step));
                
                %array containing integration steps (same sign as tk)
                ii = ones(n,1)*int_step*(tk/abs(tk));
                
                %check residual iteration step (i.e. remaining fraction of int_step)
                int_step_res = rem(tk,int_step);
                
                %adjust the total number of iterations and the array of iteration steps
                if (int_step_res ~= 0)
                    n = n + 1;
                    ii = [ii; int_step_res];
                end
                
                %numerical integration steps (i.e. re-calculation of satellite positions from toe to tk)
                pos = [X Y Z];
                vel = [Xv Yv Zv];
                acc = [Xa Ya Za];
                
                for s = 1 : n
                    
                    %Runge-Kutta numerical integration algorithm
                    %
                    %step 1
                    pos1 = pos;
                    vel1 = vel;
                    [pos1_dot, vel1_dot] = satellite_motion_diff_eq(pos1, vel1, acc, orbital_p.ELL.A, orbital_p.GM, sys_str.J2, orbital_p.OMEGAE_DOT);
                    %
                    %step 2
                    pos2 = pos + pos1_dot*ii(s)/2;
                    vel2 = vel + vel1_dot*ii(s)/2;
                    [pos2_dot, vel2_dot] = satellite_motion_diff_eq(pos2, vel2, acc, orbital_p.ELL.A, orbital_p.GM, sys_str.J2, orbital_p.OMEGAE_DOT);
                    %
                    %step 3
                    pos3 = pos + pos2_dot*ii(s)/2;
                    vel3 = vel + vel2_dot*ii(s)/2;
                    [pos3_dot, vel3_dot] = satellite_motion_diff_eq(pos3, vel3, acc, orbital_p.ELL.A, orbital_p.GM, sys_str.J2, orbital_p.OMEGAE_DOT);
                    %
                    %step 4
                    pos4 = pos + pos3_dot*ii(s);
                    vel4 = vel + vel3_dot*ii(s);
                    [pos4_dot, vel4_dot] = satellite_motion_diff_eq(pos4, vel4, acc, orbital_p.ELL.A, orbital_p.GM, sys_str.J2, orbital_p.OMEGAE_DOT);
                    %
                    %final position and velocity
                    pos = pos + (pos1_dot + 2*pos2_dot + 2*pos3_dot + pos4_dot)*ii(s)/6;
                    vel = vel + (vel1_dot + 2*vel2_dot + 2*vel3_dot + vel4_dot)*ii(s)/6;
                end
                
                %transformation from PZ-90.02 to WGS-84 (G1150)
                satp = zeros(3,1);
                satp(1,1) = pos(1) - 0.36;
                satp(2,1) = pos(2) + 0.08;
                satp(3,1) = pos(3) + 0.18;
                
                %satellite velocity
                satv = zeros(3,1);
                satv(1,1) = vel(1);
                satv(2,1) = vel(2);
                satv(3,1) = vel(3);
            end
            
        end
        
        function [az, el] = computeAzimuthElevationXS(XS, XR)
            %   Compute Azimuth and elevation of the staellite
            %
            % INPUT
            %   XS = positions of satellite [n_epoch x 1]
            %   XR = positions of reciever [n_epoch x 1] (optional, non static case)
            %
            % OUTPUT
            %   az = Azimuths of satellite [n_epoch x 1]     [deg]
            %   el = Elevations of satellite [n_epoch x 1]   [deg]
            %   at time of travel
            %
            % SYNTAX
            %   [az, el] = this.computeAzimuthElevationXS(XS)
            n_epoch = size(XS,1);
            if size(XS, 3) == 3                
                n_sat = size(XS,2);
            else
                n_sat = 1;
            end

            [az, el] = deal(zeros(n_epoch, n_sat));
            
            if n_sat > 1 % multi-satellite
                XR = reshape(XR, size(XR,1), 1, size(XR,2));
                XR = repmat(XR, iif(size(XR,1) == 1, size(XS,1), 1), size(XS, 2), 1);
            end
            
            if n_sat > 1 % multi-satellite
                [phi, lam] = cart2geod(serialize(XR(:,:,1)), serialize(XR(:,:,2)), serialize(XR(:,:,3)));
            else
                [phi, lam] = cart2geod(XR(:,1), XR(:,2), XR(:,3));
            end
            XSR = XS - XR; %%% sats orbit with origin in receiver
            XSR = reshape(XSR, n_epoch * n_sat, 3);
            
            e_unit = [-sin(lam)            cos(lam)           zeros(size(lam)) ]; % East unit vector
            n_unit = [-sin(phi).*cos(lam) -sin(phi).*sin(lam) cos(phi)]; % North unit vector
            u_unit = [ cos(phi).*cos(lam)  cos(phi).*sin(lam) sin(phi)]; % Up unit vector
            
            e = sum(e_unit .* XSR, 2);
            n = sum(n_unit .* XSR, 2);
            u = sum(u_unit .* XSR, 2);
            
            hor_dist = sqrt( e.^2 + n.^2);
            
            zero_idx = hor_dist < 1.e-20;
                
            az(zero_idx) = 0;
            el(zero_idx) = 90;
            
            az(~zero_idx) = atan2d(e(~zero_idx), n(~zero_idx));
            el(~zero_idx) = atan2d(u(~zero_idx), hor_dist(~zero_idx));
        end
        function coeff = fastLI(x_pred)
            % SYNTAX:
            %   coeff = fastLI(y_obs, x_pred);
            %
            % INPUT:
            %   y_obs  = row-vector of [1 x p_deg + 1] data values
            %   y_pred = row-vector of [1 x numel(x_pred)], where interpolation is to be found (could be a single value)
            %
            % OUTPUT:
            %   y_pred = a row-vector of interpolated y-values
            %
            % DESCRIPTION:
            %   Lagrange interpolation algorithm supposing y_obs regularly sampled
            %   The degree of the polinomial (p_deg) is equivalent to the number of
            %   element of y_obs - 1
            %
            % ORIGINAL CODE:
            %   Author: Dmitry Pelinovsky
            %   Available online at: http://dmpeli.mcmaster.ca/Matlab/Math4Q3/Lecture2-1/LagrangeInter.m
            
            n_obs = 11;
            n_pred = numel(x_pred);
            
            %y_pred = zeros(size(x_pred));
            coeff = ones(n_obs, n_pred);
            for i = 1 : n_pred
                for k = 1 : n_obs
                    for kk = 1 : (k-1) % start the inner loop through the data values for x (if k = 0 this loop is not executed)
                        coeff(kk, i) = coeff(kk, i) .* ((x_pred(i) - k) / (kk - k)); % see the Lagrange interpolating polynomials
                    end % end of the inner loop
                    
                    for kk = k+1 : n_obs % start the inner loop through the data values (if k = n this loop is not executed)
                        coeff(kk, i) = coeff(kk, i) .* ((x_pred(i) - k) / (kk - k));
                    end % end of the inner loop
                end
            end
        end
    end
    
    % ==================================================================================================================================================
    %% SHOW
    % ==================================================================================================================================================
    methods
        function fh = showOrbitsAvailability(this, start_time, stop_time, flag_no_clock)
            % Basic visualization of the orbits availability for the current
            % processing time
            %
            % SYNTAX:
            %    this.showOrbitsAvailability(<start_date>, <stop_time>)
            
            if nargin == 1
                start_time = Core.getState.getSessionsStartExt;
                stop_time = Core.getState.getSessionsStopExt;
            end
            
            if nargin == 2
                stop_time = start_time.last();
                start_time = start_time.first();
            end
            
            if nargin <= 3 || isempty(flag_no_clock)
                flag_no_clock = false;
            end
            
            this.initSession(start_time, stop_time);
            
            fh = figure;
            subplot(16,1,1:6); Core_UI.addBeautifyMenu(fh); drawnow
            if not(isempty(this.getClockTime))
                t_clock = this.getClockTime.getMatlabTime;
                imagesc(t_clock, 1 : size(this.clock, 2), not(isnan(zero2nan(this.clock)))');
            end
            xlim([start_time.getMatlabTime, stop_time.getMatlabTime]);
            setTimeTicks();
            ax(1) = gca;
            title(sprintf('Clock availability\\fontsize{5} \n'));
            ylabel('Satellites');

            subplot(16,1,9:16);
            cmap = [0.2 0.2 0.2; 1 0.5 0.1; 0.3 1 0.3];
            coord_validity = not(isnan(zero2nan(this.coord(:,:,1))));
            if not(isempty(coord_validity))
                coord_validity = coord_validity + flagShrink(coord_validity, 5);
                t_clock = this.getCoordTime.getMatlabTime;
                imagesc(t_clock, 1 : size(this.coord, 2), coord_validity');
            end
            xlim([start_time.getMatlabTime, stop_time.getMatlabTime]);
            colormap(cmap);
            setTimeTicks();
            ax(2) = gca;
            title(sprintf('Coordinates availability\\fontsize{5} \n'));
            ylabel('Satellites');
            
            caxis([0 2]);
            cb = colorbar('Location', 'SouthOutside');
            cb.Ticks  = (1:2:5)/3;
            cb.TickLabels = {'No data', 'Polynomial border', 'Good data'};

            linkaxes(ax);
            xlim([start_time.getMatlabTime, stop_time.getMatlabTime]);
            
            fig_name = sprintf('Orbits_availability');
            fh.UserData = struct('fig_name', fig_name);
            Core_UI.beautifyFig(fh);
            Core_UI.addExportMenu(fh);
            Core_UI.addBeautifyMenu(fh);
        end
    end
end
