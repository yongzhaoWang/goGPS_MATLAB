%   CLASS GNSS_Station
% =========================================================================
%
%
%   Class to store receiver data (observations, and characteristics)
%
% EXAMPLE
%   trg = GNSS_Station();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Receiver

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b8
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti, Giulio Tagliaferro ...
%  Contributors:
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
classdef GNSS_Station < handle
    properties (SetAccess = private, GetAccess = private)
        creation_time = GPS_Time(now); % object creation time
    end

    properties (SetAccess = public, GetAccess = public)
        marker_name    % marker name
        marker_type    % marker type
        number         % receiver number
        type           % receiver type
        version        % receiver version
        observer       % name of observer
        agency         % name of agency

        % ANTENNA ----------------------------------

        ant_serial     % antenna number
        ant_type       % antenna type
        ant_delta_h    % antenna height from the ground [m]
        ant_delta_en   % antenna east/north offset from the ground [m]
        ant_mp         % antenna multipath mitigation model

        static         % static or dynamic receiver 1: static 0: dynamic

        work           % handle to receiver Work Space
        old_work       % handle to the old Work Space used for repair CS (contains just few epochs)
        out            % handle to receiver outputs

        cc             % constallation collector contains information on the Constellation used
        state          % handle to the state object

        log            % handle to the log object
        w_bar          % handle to the wait bar object
    end

    % ==================================================================================================================================================
    %% PROPERTIES PLOTS
    % ==================================================================================================================================================

    properties
        slant_filter_win = 0; % used in some visualization represente the knot distance of splines used for filtering
    end

    % ==================================================================================================================================================
    %% METHODS INIT - CLEAN - RESET - REM - IMPORT
    % ==================================================================================================================================================
    methods
        function this = GNSS_Station(flag_static)
            % Creator method
            %
            % INPUT
            %   flag_static  flag is static [ boolean ]
            %
            % SYNTAX
            %   this = GNSS_Static(static)
            this.work = Receiver_Work_Space(this);
            this.out = Receiver_Output(this);
            if nargin >= 2 && ~isempty(flag_static)
                this.static = logical(flag_static);
            end
            this.init();
            this.resetInfo();
        end

        function importRinexLegacy(this, rinex_file_name, rate, sys_c_list, otype_list)
            % Select the files to be imported
            %
            % INPUT
            %   rinex_file_name     path to a RINEX file
            %   rate                import rate
            %   sys_c_list          list of char specifing the satellite system to read
            %
            % SYNTAX
            %   this.importRinexLegacy(rinex_file_name, rate)
           if ~isempty(rinex_file_name) && (exist(rinex_file_name, 'file') == 2)
                this.work.rinex_file_name = rinex_file_name;
            else
                this.work.rinex_file_name = '';
            end
            this.work.load(rate, sys_c_list, otype_list);
            this.work.out_start_time = this.work.time.first;
            this.work.out_stop_time = this.work.time.last;
        end

        function importRinexes(this, rin_list, time_start, time_stop, rate, sys_c_list, otype_list)
            % Select the files to be imported
            %
            % INPUT
            %   rin_list      object containing the list of rinex to load [ File_Rinex ]
            %   time_start    first epoch to load [ GPS_Time ]
            %   time_stop     last epoch to load [ GPS_Time ]
            %   rate          import rate [s]
            %   sys_c_list          list of char specifing the satellite system to read
            %
            % SYNTAX
            %   this.importRinexes(rin_list, time_start, time_stop, rate)
            this.work.importRinexFileList(rin_list, time_start, time_stop, rate, sys_c_list, otype_list);
        end
        
        function synthetizeWork(this, xyz, time, system)
            % synththize observations
            %
            % SYNTAX:
            %   this.synthetizeWork(xyz)
            work = this.work;
            work.xyz = xyz;
            work.time = time;
            
            % two pseudoranges and two phases per satellites
            cc = Core.getConstellationCollector;
            work.go_id  = repmat(cc.index,4,1);
            work.prn    = repmat(cc.prn,4,1);
            work.system = repmat(cc.system,1,4);
            
            o_code1 = [];
            o_code2 = [];
            wl1 = [];
            wl2 = [];
            for ss_ch = unique(work.system)
                n_sat = sum(cc.system == ss_ch);
                sys = cc.getSys(ss_ch);
                if ss_ch == 'G'
                    trk = 2;
                else
                    trk = 1;
                end
                o_code1 = [o_code1; repmat([sys.CODE_RIN3_2BAND(1) sys.CODE_RIN3_ATTRIB{1}(trk)],n_sat,1)];
                o_code2 = [o_code2; repmat([sys.CODE_RIN3_2BAND(2) sys.CODE_RIN3_ATTRIB{2}(trk)],n_sat,1)];
                if ss_ch ~= 'R'
                    wl1 = [wl1; repmat(sys.L_VEC(1),n_sat,1)];
                    wl2 = [wl2; repmat(sys.L_VEC(2),n_sat,1)];
                else
                    wl1 = [wl1; sys.L_VEC(ss.PRN2IDCH(min(sys.PRN, sys.N_SAT))')];
                    wl2 = [wl2; sys.L_VEC(size(sys.L_VEC, 1) + ss.PRN2IDCH(min(sys.PRN, sys.N_SAT))')];
                end
            end
            
            work.obs_code = [[repmat('C',size(o_code1,1),1) o_code1]; [repmat('L',size(o_code1,1),1) o_code1]; [repmat('C',size(o_code2,1),1) o_code2]; [repmat('L',size(o_code2,1),1) o_code2]];
            work.wl = [wl1;wl1;wl2;wl2];
            work.f_id = [ones(size(wl1));ones(size(wl1));2*ones(size(wl2));2*ones(size(wl2))];
            work.active_ids = true(size(work.wl));
            work.obs = zeros(size(work.obs_code,1),work.time.length);
            work.sat.avail_index = true(size(work.obs'));
            work.updateAzimuthElevation();
            work.sat.avail_index = work.sat.el > 2;
            work.n_spe = sum(work.sat.avail_index,2);
            work.updateAzimuthElevation();
            work.updateErrIono();
            err_iono = work.sat.err_iono; % consider add random fluctuataions
            work.updateErrTropo();
            L1SQ = GPS_SS.L_VEC(1)^2;
            rec_set = Receiver_Settings;
            dt = cumsum(randn(1,work.time.length));
            for i = 1:4
                [range] = work.getSyntObs( unique(work.go_id));
                work.sat.tot = nan2zero((zero2nan(range)'- repmat(dt',1,size(range,1),1) + err_iono*L1SQ)) / Core_Utils.V_LIGHT;
            end

            for ss_ch = unique(work.system)
                sys = cc.getSys(ss_ch);
                idx_sat = cc.index(cc.system == ss_ch);
                [range] = work.getSyntObs( unique(work.go_id(work.system == ss_ch)));
                % systethize C1
                idx_obs = work.obs_code(:,1) == 'C' & work.system' == ss_ch & work.obs_code(:,2) == sys.CODE_RIN3_2BAND(1);
                noise = rec_set.getStd(ss_ch, work.obs_code(find(idx_obs,1,'first'),:));
                work.obs(idx_obs,:) = nan2zero(zero2nan(range) - zero2nan(err_iono(:,idx_sat)'.*repmat(work.wl(idx_obs),1,work.time.length).^2) + rand(1) + noise * randn(size(range))); 
                % systethize C2
                idx_obs = work.obs_code(:,1) == 'C' & work.system' == ss_ch & work.obs_code(:,2) == sys.CODE_RIN3_2BAND(2);
                noise = rec_set.getStd(ss_ch, work.obs_code(find(idx_obs,1,'first'),:));
                work.obs(idx_obs,:) = nan2zero(zero2nan(range) - zero2nan(err_iono(:,idx_sat)'.*repmat(work.wl(idx_obs),1,work.time.length).^2) + rand(1) + noise * randn(size(range))); 
                % systethize L1
                idx_obs = work.obs_code(:,1) == 'L' & work.system' == ss_ch & work.obs_code(:,2) == sys.CODE_RIN3_2BAND(1);
                noise = rec_set.getStd(ss_ch, work.obs_code(find(idx_obs,1,'first'),:));
                work.obs(idx_obs,:) = nan2zero(zero2nan(range) + zero2nan(err_iono(:,idx_sat)'.*repmat(work.wl(idx_obs),1,work.time.length).^2) + rand(1) + noise * randn(size(range))); 
                
                % systethize L2
                idx_obs = work.obs_code(:,1) == 'L' & work.system' == ss_ch & work.obs_code(:,2) == sys.CODE_RIN3_2BAND(2);
                noise = rec_set.getStd(ss_ch, work.obs_code(find(idx_obs,1,'first'),:));
                work.obs(idx_obs,:) = nan2zero(zero2nan(range) + zero2nan(err_iono(:,idx_sat)'.*repmat(work.wl(idx_obs),1,work.time.length).^2) + rand(1) + noise * randn(size(range)));
            end
            work.group_delay_status= true;
            work.dts_delay_status= true;
            work.sh_delay_status= true;
            work.pcv_delay_status= true;
            work.ol_delay_status= true;
            work.pt_delay_status= true;
            work.pw_delay_status= true;
            work.et_delay_status= true;
            work.hoi_delay_status= true;
            work.atm_load_delay_status= true;
            % add dt
            work.obs = nan2zero(zero2nan(work.obs) + repmat(dt,size(work.obs,1),1));
            idx_phase = work.obs_code(:,1) == 'L';
            work.obs(idx_phase,:) = work.obs(idx_phase,:)./repmat(work.wl(idx_phase),1,size( work.obs(idx_phase,:),2));
            

        end

        function clearHandles(this)
            % Clear handles
            %
            % SYNTAX
            %   this.clearHandles();

            this.log = [];

            this.w_bar = [];
        end
        
        function initHandles(this)
            % Reload handles
            % 
            % SYNTAX
            %   this.initHandles
            
            this.log = Core.getLogger();

            this.w_bar = Go_Wait_Bar.getInstance();
        end
        
        function init(this)
            % Reset handles
            %
            % SYNTAX
            %   this.init();

            this.log = Core.getLogger();

            this.w_bar = Go_Wait_Bar.getInstance();
            this.work = Receiver_Work_Space(this);
        end

        function resetInfo(this)
            % Reset information about receiver (name, type, number...)
            %
            % SYNTAX
            %   this.reset()
            this.marker_name  = 'unknown';  % marker name
            this.marker_type  = '';       % marker type
            this.number   = '000';
            this.type     = 'unknown';
            this.version  = '000';
        end

        function resetWork(sta_list)
            % Reset handle to work object
            %
            % SYNTAX
            %   this.resetWork()
            for r = 1 : numel(sta_list)
                sta_list(r).work.resetWorkSpace();
                sta_list(r).old_work = Receiver_Work_Space(sta_list(r));
            end
        end

        function resetOut(sta_list)
            % Reset handle to output object
            %
            % SYNTAX
            %   this.resetOut()
            for r = 1 : numel(sta_list)
                sta_list(r).out = Receiver_Output(sta_list(r));
            end
        end

        function netPrePro(sta_list)
            % EXPERIMENTAL pre processing multi-receiver
            % Perform multi-receiver pre-processing
            %  - outlier detection,
            %  - realignment of "old" phases for ambiguity passing + repair
            %
            % INPUT
            %   sta_list    list of receiver
            %
            % SYNTAX
            %   sta_list.netProPro()

            flag_ct_apply = true; % Apply the estimated clock to the data and repeat outlier detection
            realign_ph = true;
            out_det = true;

            show_fig = false;
            % Prepare data in a unique structure

            work_list = [sta_list(~sta_list.isEmptyWork_mr).work];
            if numel(work_list) > 1 && (show_fig || out_det)

                [t_sync, id_rsync] = Receiver_Commons.getSyncTimeExpanded(work_list);
                %id_rsync(any(isnan(zero2nan(id_rsync)')), :) = [];

                n_epochs = size(id_rsync, 1);
                n_rec = numel(work_list);

                clear dt_red ph_red id_ph_red
                for r = 1 : n_rec
                    %work_list(r).keepBestTracking();
                    [dt_red{r}, ph_red{r}, id_ph_red{r}] = work_list(r).getReducedPhases();
                    dt_red{r} = dt_red{r}(noNaN(id_rsync(:, r)), :);
                    %ph_red{r} = bsxfun(@rdivide, ph_red{r}(noNaN(id_rsync(:, r)), :), work_list(r).wl(id_ph_red{r})');
                end

                % Get all SS present in the receivers
                all_ss = unique([work_list.system]);
                for sys_c = all_ss
                    n_obs = 0;
                    prn_list = [];
                    bands = '';
                    % Each phase will be added
                    id_sat = {};
                    for r = 1 : n_rec
                        obs_code = work_list(r).getAvailableObsCode('L', sys_c);
                        bands = [bands; obs_code(:,2:3)];
                        n_obs = n_obs + size(obs_code, 1);
                        id_ph{r} = work_list(r).findObservableByFlag('L');
                        id_sat{r} = work_list(r).system(id_ph{r}) == sys_c;
                        prn_list = unique([prn_list; work_list(r).prn(id_ph{r}(id_sat{r}))]);
                    end
                    n_sat = numel(prn_list);
                    tracking = unique(bands(:,2));
                    bands = unique(bands(:,1));
                    n_bands = numel(bands);

                    all_ph_red = nan(n_epochs, n_sat * n_bands, n_rec);
                    all_dph_red = nan(n_rec, n_epochs, n_sat * n_bands);
                    for r = 1 : n_rec
                        % Find all the phases
                        id = work_list(r).findObservableByFlag('L', sys_c);
                        [id_ok, ~, id_red] = intersect(id, id_ph_red{r});
                        for i = 1 : numel(id_ok)
                            [~, ~, p_list] = intersect(work_list(r).prn(id_ok(i)), prn_list);
                            [~, ~, b_list] = intersect(work_list(r).obs_code(id_ok(i), 2), bands);
                            
                            sid = repmat(p_list, n_bands, 1) + serialize(repmat(numel(prn_list) * (b_list - 1)', numel(p_list), 1)); % < -this only works if all band are available on all satellites, otherwise it will crash, to be fixed
                            all_ph_red(not(isnan(id_rsync(:,r))), sid, r) = zero2nan(ph_red{r}(:, id_red(i)));
                            tmp = Core_Utils.diffAndPred(all_ph_red(:, sid, r));
                            tmp = bsxfun(@minus, tmp, strongMean(tmp,0.95, 0.95, 2));
                            id_ko = work_list(r).sat.outliers_ph_by_ph(noNaN(id_rsync(:, r)), id_sat{r}(id_red(i))) | work_list(r).sat.cycle_slip_ph_by_ph(noNaN(id_rsync(:, r)), id_sat{r}(id_red(i)));
                            tmp(not(isnan(id_rsync(:,r)))) = tmp(not(isnan(id_rsync(:,r)))) .* (1-id_ko);
                            all_dph_red(r, :, sid) = zero2nan(tmp);
                        end
                    end

                    if show_fig || out_det
                        % Estimate common term from data
                        dct = nan(size(all_dph_red,2), size(all_dph_red,3));
                        for s = 1 : size(dct, 2)
                            dct(:,s) = strongMean(all_dph_red(:,:,s));
                        end
                        dct = squeeze(median(all_dph_red, 1, 'omitnan'));
                        dct(sum(~isnan(zero2nan(all_dph_red))) <= 1) = 0;

                        id_even = squeeze(sum(~isnan(zero2nan(all_dph_red))) == 2);
                        if any(id_even(:))
                            % find the observation that is closer to zero
                            [tmp, id] = min(abs(zero2nan(all_dph_red)));
                            tmp_min = all_dph_red(id(:) + 2*(0 : (numel(dct) -1))');

                            % when I have 2 observations choose the observation closer to zero
                            dct(id_even) = tmp_min(id_even);
                        end
                        dct = nan2zero(dct);
                        ct = cumsum(dct);
                        ct(dct == 0) = nan;
                    end

                    if out_det
                        if flag_ct_apply
                            [ph, wl] = work_list(r).getPhases(sys_c);
                            cur_id_ph = id_ph{r}(id_sat{r});
                            
                            [id_ok, ~, id_red] = intersect(cur_id_ph, id_ph_red{r});
                            for i = 1 : numel(id_ok)
                                [~, ~, p_list] = intersect(work_list(r).prn(id_ok(i)), prn_list);
                                [~, ~, b_list] = intersect(work_list(r).obs_code(id_ok(i), 2), bands);
                                
                                sid = repmat(p_list, n_bands, 1) + serialize(repmat(numel(prn_list) * (b_list - 1)', numel(p_list), 1)); % < -this only works if all band are available on all satellites, otherwise it will crash, to be fixed
                                ph(:, i) = ph(:, i) + ct(noNaN(id_rsync(:,r)), sid);
                            end
                            work_list(r).setPhases(ph, wl, cur_id_ph)
                        end
                        
                        for r = 1 : n_rec
                            id = work_list(r).findObservableByFlag('L', sys_c);
                            [id_ok, ~, id_red] = intersect(id, id_ph_red{r});
                            id_ko = false(size(work_list(r).sat.outliers_ph_by_ph));
                            for i = 1 : numel(id_ok)
                                [~, ~, p_list] = intersect(work_list(r).prn(id_ok(i)), prn_list);
                                [~, ~, b_list] = intersect(work_list(r).obs_code(id_ok(i), 2), bands);
                                
                                sid = repmat(p_list, n_bands, 1) + serialize(repmat(numel(prn_list) * (b_list - 1)', numel(p_list), 1)); % < -this only works if all band are available on all satellites, otherwise it will crash, to be fixed
                                
                                sensor = abs((squeeze(all_dph_red(r,not(isnan(id_rsync(:,r))),sid))' - zero2nan(dct(not(isnan(id_rsync(:,r))), sid))) .* (abs(dct(not(isnan(id_rsync(:,r))), sid)) > 0)) > 0.05;
                                
                                % V0
                                % id_ph = work_list(r).findObservableByFlag('L', sys_c);
                                % for b = b_list'
                                %     for p = 1: numel(p_list)
                                %         sid = p + numel(prn_list) * (b - 1);
                                %
                                %         id_obs = work_list(r).findObservableByFlag(['L' bands(b)], sys_c, prn_list(p_list(p)));
                                %         [~, ~, ido] = intersect(id_obs, id_ph);
                                %         work_list(r).sat.outliers_ph_by_ph(id_rsync(:,r), ido) = work_list(r).sat.outliers_ph_by_ph(id_rsync(:,r), ido) | ;
                                %     end
                                % end
                                
                                % V1 (improve V0)
                                id_ko(noNaN(id_rsync(:,r)),id_sat{r}(id_red(i))) = sensor;
                            end
                            work_list(r).addOutliers(id_ko, true);
                        end
                    end
                end

                % Realign phases with the past
                % (useful when work_list(r).state.flag_amb_pass)
                if realign_ph
                    state = Core.getState;
                    for r = 1 : n_rec
                        if state.flag_amb_pass && ~isempty(work_list(r).parent.old_work) && ~work_list(r).parent.old_work.isEmpty
                            t_new = round(work_list(r).parent.work.time.getRefTime(work_list(r).parent.old_work.time.first.getMatlabTime) * 1e7) / 1e7;
                            t_old = round(work_list(r).parent.old_work.time.getRefTime(work_list(r).parent.old_work.time.first.getMatlabTime) * 1e7) / 1e7;
                            [~, id_new, id_old] = intersect(t_new, t_old);

                            if ~isempty(id_new)
                                [ph, wl, lid_ph] = work_list(r).getPhases();
                                tmp = ph - work_list(r).getSyntPhases;
                                id_ph = find(lid_ph);
                                for i = 1 : length(id_ph)
                                    [amb_off, old_ph, old_synt] = work_list(r).parent.old_work.getLastRepair(work_list(r).go_id(id_ph(i)), work_list(r).obs_code(id_ph(i),2:3));
                                    if ~isempty(old_ph)
                                        ph_diff = (tmp(id_new, i) / wl(i) - (old_ph(id_old) - old_synt(id_old)));
                                        amb_off_emp = median(round(ph_diff / state.getCycleSlipThr()), 'omitnan') * state.getCycleSlipThr;
                                        if ~isempty(amb_off_emp) && amb_off_emp ~= 0
                                            ph(:,i) = ph(:,i) - amb_off_emp * wl(i);
                                        end
                                    end
                                end
                                work_list(r).setPhases(ph, wl, lid_ph);
                            end
                        end
                    end
                end

                % plots
                if show_fig
                    for sys_c = all_ss

                        for r = 1 : n_rec
                            figure; plot(squeeze(all_dph_red(r,:,:)) - dct); title(sprintf('Receiver %d diff', r));
                            figure; clf;
                            [tmp, tmp_trend, tmp_jmp] = work_list(r).flattenPhases(squeeze(all_ph_red(:, :, r)) - cumsum(dct));
                            plot(all_ph_red(:, :, r) - tmp_trend + tmp_jmp - repmat(strongMean(all_ph_red(:, :, r) - tmp_trend + tmp_jmp), size(tmp_trend, 1), 1));
                            title(sprintf('Receiver (Possible Repair) %d full', r));
                            dockAllFigures;
                            figure; clf;
                            plot(all_ph_red(:, :, r) - tmp_trend - repmat(strongMean(all_ph_red(:, :, r) - tmp_trend), size(tmp_trend, 1), 1));
                            title(sprintf('Receiver %d full', r));
                            dockAllFigures;
                            try
                                ww = work_list(r).parent.old_work; ph_red = ww.getPhases() - ww.getSyntPhases(); figure; plot(ww.time.getMatlabTime, ph_red ./ ww.wl(1), '.-k')
                                ww = work_list(r); ph_red = ww.getPhases() - ww.getSyntPhases(); hold on; plot(ww.time.getMatlabTime, ph_red ./ ww.wl(1))
                            catch
                            end
                            dockAllFigures;
                        end
                    end
                end
            end
        end
    end
    
    %% METHODS ADVANCED
    % ==================================================================================================================================================
    methods
        function updateMultiPath(sta_list, day_span, mp_mode, sys_grp)
            sta_list = sta_list(~sta_list.isEmpty);
            log = Core.getLogger();
            if (nargin < 4) || isempty(sys_grp)
                sys_grp = struct('G', 'G', 'R', 'R', 'E', 'E', 'J', 'J', 'C', 'C', 'I', 'I', 'S', 'S');
            end
            for rec = sta_list(:)'
                log.addMarkedMessage(sprintf('Updating multipath corrections for "%s"', rec.getMarkerName4Ch));
                % If resduals ph by ph are available on out use them
                if nargin >= 2 && ~isempty(day_span)
                    if isempty(rec(1).out.sat.res)
                        ant_mp = rec.work.computeMultiPath(sys_grp, [], day_span);
                    else
                        ant_mp = rec.out.computeMultiPath(sys_grp, [], day_span);
                    end
                else
                    if isempty(rec(1).out.sat.res)
                        ant_mp = rec.work.computeMultiPath(sys_grp);
                    else
                        ant_mp = rec.out.computeMultiPath(sys_grp);
                    end
                end
                
                % rec.ant_mp = rec.ant_mp + ant_mp;
                flag_update = false;
                if isempty(rec.ant_mp) && (nargin >= 3 && not(isempty(mp_mode)))
                    rec.ant_mp = rec.getAntennaMultiPath();
                end
                if (isempty(rec.work) || isempty(rec.work.ant_mp)) && ~(nargin >= 3 && not(isempty(mp_mode)))
                    % Zerniche multipath is not yet in the receiver
                    rec.ant_mp = ant_mp;
                    flag_update = true;
                else
                    applied_ant = rec.work.ant_mp;
                    % Zerniche multipath is already in the receiver
                    ant_mp_bk = rec.ant_mp;
                    try
                        % Get the satellite systems available in the zerniche multipath struct
                        sys_c_list = cell2mat(fields(ant_mp)');
                        sys_c_list = intersect(sys_c_list, 'GREJCIS');
                        for sys_c = sys_c_list
                            trk_list = fields(ant_mp.(sys_c))';
                            if ~isfield(rec.ant_mp, sys_c)
                                % This constellation is not present into the old Zernike MultiPath set of coefficients
                                rec.ant_mp.(sys_c) = ant_mp.(sys_c);
                            else
                                for trk = trk_list
                                    trk = strrep(trk{1}, ' ', '_'); % Spaces are not supported in structures
                                    short_trk = trk;
                                    if numel(short_trk) > 3 && all(short_trk(4:end) == '_')
                                        short_trk = short_trk(1:3);
                                    end
                                    if (~isfield(rec.ant_mp.(sys_c), trk) && ~isfield(rec.ant_mp.(sys_c), short_trk)) || ((trk(end) == 'I') && Core.getCurrentSettings.flag_rec_mp == 0)
                                        % This tracking frequency is not present into the old Zernike MultiPath set of coefficients
                                        rec.ant_mp.(sys_c).(short_trk) = ant_mp.(sys_c).(trk);
                                    else
                                        % Get the old map that was applied
                                        if nargin >= 3 && not(isempty(mp_mode))
                                            applied_ant = rec.ant_mp;
                                            if isempty(applied_ant)
                                                applied_ant = rec.getAntennaMultiPath();
                                            end
                                            if isfield(applied_ant.(sys_c), trk)
                                                app_trk = trk; % longname
                                            else
                                                app_trk = short_trk; % shortname
                                            end
                                            switch mp_mode
                                                case 0, applied_map = 0;
                                                case 1, applied_map = applied_ant.(sys_c).(app_trk).z_map;
                                                case 2, applied_map = applied_ant.(sys_c).(app_trk).r_map;
                                                case 3, applied_map = applied_ant.(sys_c).(app_trk).g_map;
                                                case 4, applied_map = applied_ant.(sys_c).(app_trk).c_map;
                                                case 5, applied_map = applied_ant.(sys_c).(app_trk).g1_map;
                                                case 6, applied_map = applied_ant.(sys_c).(app_trk).c1_map;
                                            end
                                        else
                                            if ~isfield(applied_ant.(sys_c), trk)
                                                app_trk = trk; % longname
                                            else
                                                app_trk = short_trk; % shortname
                                            end
                                            switch Core.getCurrentSettings.flag_rec_mp
                                                case 0, applied_map = 0;
                                                case 1, applied_map = applied_ant.(sys_c).(app_trk).z_map;
                                                case 2, applied_map = applied_ant.(sys_c).(app_trk).r_map;
                                                case 3, applied_map = applied_ant.(sys_c).(app_trk).g_map;
                                                case 4, applied_map = applied_ant.(sys_c).(app_trk).c_map;
                                                case 5, applied_map = applied_ant.(sys_c).(app_trk).g1_map;
                                                case 6, applied_map = applied_ant.(sys_c).(app_trk).c1_map;
                                            end
                                        end
                                        if isfield(rec.ant_mp.(sys_c), trk) && numel(trk) > 3
                                            % I prefer short names
                                            rec.ant_mp.(sys_c).(short_trk) = rec.ant_mp.(sys_c).(trk); % copy short name
                                            rec.ant_mp.(sys_c) = rmfield(rec.ant_mp.(sys_c), trk); % remove long name;
                                        end
                                        rec.ant_mp.(sys_c).(short_trk).z_map  = Core_Utils.resize2(applied_map, size(ant_mp.(sys_c).(trk).z_map)) + ant_mp.(sys_c).(trk).z_map;
                                        rec.ant_mp.(sys_c).(short_trk).r_map  = Core_Utils.resize2(applied_map, size(ant_mp.(sys_c).(trk).r_map)) + ant_mp.(sys_c).(trk).r_map;
                                        rec.ant_mp.(sys_c).(short_trk).g_map  = Core_Utils.resize2(applied_map, size(ant_mp.(sys_c).(trk).g_map)) + ant_mp.(sys_c).(trk).g_map;
                                        rec.ant_mp.(sys_c).(short_trk).c_map  = Core_Utils.resize2(applied_map, size(ant_mp.(sys_c).(trk).c_map)) + ant_mp.(sys_c).(trk).c_map;
                                        rec.ant_mp.(sys_c).(short_trk).g1_map = Core_Utils.resize2(applied_map, size(ant_mp.(sys_c).(trk).g1_map)) + ant_mp.(sys_c).(trk).g1_map;
                                        rec.ant_mp.(sys_c).(short_trk).c1_map = Core_Utils.resize2(applied_map, size(ant_mp.(sys_c).(trk).c1_map)) + ant_mp.(sys_c).(trk).c1_map;
                                    end
                                end
                            end
                        end
                        rec.ant_mp.time_lim = ant_mp.time_lim;
                        flag_update = true;
                    catch ex
                        %Core_Utils.printEx(ex);
                        log.addError(sprintf('The new multipath model for "%s" is not compatible with the previous one\nDelete or do not load the old maps', rec.getMarkerName4Ch));
                        % if any error arises this set is not compatible with the previous one
                        % e.g. it could have different maximum degree, or different frequencies
                        rec.ant_mp = ant_mp_bk;
                    end
                end
                if flag_update && Core.getCurrentSettings.flag_rec_mp > 0
                    % If I add the zerniche polynomials to the receiver
                    % I also correct the observations with the new set
                    if isempty(rec(1).out.sat.res_ph_by_ph) % if I work on out do not apply MP
                        rec.work.applyMP();
                    end
                    log.addStatusOk(sprintf('Update completed for "%s"', rec.getMarkerName4Ch));
                end
            end
        end
        
        function zhd_corr = getZenithDelayCorrection(this, time, coo, pressure)
            coo_gnss = this.getPos.getMedianPos();
            [pr_gnss, ~, ~, p_time, id_sync] = this.getPTH_mr();            
                        
            p_time = p_time.getEpoch(id_sync).getMatlabTime;
            pr_gnss = pr_gnss(id_sync);           
            
            pressure_gnss = interp1q(p_time, pr_gnss, time.getMatlabTime);
            
            if nargin == 3
                zhd_corr = Atmosphere.getZenithDelayCorrection(coo_gnss, pressure_gnss, coo);
            else
                zhd_corr = Atmosphere.getZenithDelayCorrection(coo_gnss, pressure_gnss, coo, pressure);
            end
        end
    end
    % ==================================================================================================================================================    
    %% METHODS EXPORT
    % ==================================================================================================================================================
    
    methods
        function exportCRD(sta_list, mode, flag)
            % Export coordinates as a CRD file
            %
            % INPUT
            %   mode    it could be 'out' (dafault), 'work', ...(future modes)
            %
            % SYNTAX
            %  sta_list.fixPos(mode)
            
            %create a new RF
            rf = Core_Reference_Frame;
            if nargin < 3 || isempty(flag)
                flag = 3; % use this coordinate for prepro
            end
            if nargin <2 || isempty(mode)
                mode = 'out'; % use median out coordinates
            end
            
            for s = 1 : length(sta_list)
                if strcmpi(mode,'out') || isempty(sta_list(s).out) || (sta_list(s).out.isEmpty)
                    xyz = sta_list(s).out.getPosXYZ();
                    xyz = xyz(end,:);
                    rf.setCoo(upper(sta_list(s).getMarkerName4Ch), xyz, flag, [0 0 0], [0 0], GPS_Time([1970 1 1 0 0 0]), GPS_Time([2099 1 1 0 0 0]));
                else %if strcmpi(mode,'work') % get from work
                    xyz = sta_list(s).work.getMedianPosXYZ();
                    rf.setCoo(upper(sta_list(s).getMarkerName4Ch), xyz, flag, [0 0 0], [0 0], GPS_Time([1970 1 1 0 0 0]), GPS_Time([2099 1 1 0 0 0]));
                end
            end
            
            out_dir = Core.getState.getOutDir();
            out_file_name = fullfile(out_dir, sprintf('coordinates_%s.crd',GPS_Time.now.toString('yyyymmdd_HHMMSS')));
            log = Core.getLogger;
            log.addMarkedMessage(sprintf('Exporting coordinates to %s',out_file_name));
            rf.export(out_file_name);
            log.addStatusOk('Export completed successfully');
        end
        
        function exportAppendedCoo(sta_list, mode, out_file_prefix)
            % Export the current value of the coordinate to a text coordinate file
            % One file per receiver, containing: time_stamp: 
            %   X, Y, Z, 
            %   Cxx, Cyy, Czz, Cxy, Cxz, Cyz
            %   n_epochs, n_obs
            %
            % INPUT:
            %   mode             it could be 'out' (dafault), 'work', ...(future modes)
            %   out_file_prefix  prefix of the out file
            %
            % SYNTAX:
            %   this.exportAppendedCoo(this, <mode = 'out'>, <out_file_name>)
            
            % Remove empty receivers
            if nargin < 2 || isempty(mode)
                mode = 'out'; % use median out coordinates
            end
            flag_out = strcmpi(mode, 'out');
            
            sta_list = sta_list(~sta_list.isEmpty_mr);
            if ~isempty(sta_list)
                log = Core.getLogger;
                for rec = sta_list(:)'
                    if flag_out && not(rec.isEmptyOut_mr)
                        coo = rec.out.getPos;
                    else
                        coo = rec.work.getPos;
                    end
                     
                    if nargin < 3 || isempty(out_file_prefix)
                        out_file_name = coo.getCooOutPath();
                    else
                        out_file_name = coo.getCooOutPath(out_file_prefix);
                    end
                    coo.exportAsCoo(out_file_name)
                end
            end
        end
        
        function exportAppendedBernyCoo(sta_list, mode, out_file_prefix)
            % Export the current value of the coordinate to a CRD and OUT file as Bernese software does
            %
            % INPUT:
            %   mode             it could be 'out' (dafault), 'work', ...(future modes)
            %   out_file_prefix  prefix of the out file
            %
            % SYNTAX:
            %   this.exportAppendedBernyCoo(this, <mode = 'out'>, <out_file_name>)
            
            % Remove empty receivers
            if nargin < 2 || isempty(mode)
                mode = 'out'; % use median out coordinates
            end
            flag_out = strcmpi(mode, 'out');
            
            sta_list = sta_list(~sta_list.isEmpty_mr);
            if ~isempty(sta_list)
                log = Core.getLogger;
                for rec = sta_list(:)'
                    if flag_out && not(rec.isEmptyOut_mr)
                        coo = rec.out.getPos;
                    else
                        coo = rec.work.getPos;
                    end
                    
                    if nargin < 3 || isempty(out_file_prefix)
                        out_file_name = coo.getOutPath();
                    else
                        out_file_name = coo.getOutPath(out_file_prefix);
                    end
                    coo.exportAsBerny(out_file_name)
                end
            end
        end
        
        function exportPlainCoord(sta_list, type, out_file_name)
            % Export the current value of the coordinate to a text coordinate file
            %
            % INPUT:
            %   type    it can be ( "XYZ" | "ENU" | "Geodetic")
            %
            % SYNTAX:
            % this.exportXYZ(this, type, <out_file_name>)
            
            % Remove empty receivers
            if nargin < 2 || isempty(type)
                type = 'XYZ';
            end
            
            if upper(type(1)) == 'X' % XYZ
               type = 'XYZ';
            elseif upper(type(1)) == 'E' % ENU
               type = 'ENU';
            else % if upper(type(1)) == 'G' % Geodetic
               type = 'Geodetic';
            end
            
            sta_list = sta_list(~sta_list.isEmptyOut_mr);
            if ~isempty(sta_list)
                out_list = [sta_list.out];
                
                state = Core.getState();
                now_time = GPS_Time.now();
                if nargin < 3 || isempty(out_file_name)
                    out_dir = state.getOutDir();
                    if numel(sta_list) == 1
                        prefix = strrep(sta_list.getMarkerName, ' ', '_');
                        prefix = sprintf('%s_%s_%s', prefix, out_list.time_pos.first.toString('yyyymmdd_HHMMSS'), out_list.time_pos.last.toString('yyyymmdd_HHMMSS'));
                    else
                        prefix = strrep(state.getPrjName, ' ', '_');
                    end
                    out_file_name = fullfile(out_dir, sprintf('Pos_%s_%s_%s.txt', type, prefix, now_time.toString('yyyymmdd_HHMMSS')));
                else
                    % Add the folder if not present
                    if sum(out_file_name == filesep) == 0
                        out_dir = state.getOutDir();
                        out_file_name = fullfile(out_dir, out_file_name);
                    end
                end
                
                n_sta = numel(sta_list);
                
                % Find first and last stored epoch for each receiver
                st_time = zeros(n_sta, 1);
                en_time = zeros(n_sta, 1);
                for r = 1 : numel(out_list)
                    st_time(r) = out_list(r).time_pos.first.getMatlabTime;
                    en_time(r) = out_list(r).time_pos.last.getMatlabTime;
                end
                
                Core.getLogger.addMarkedMessage(sprintf('Exporting %s coordinates to %s', type, out_file_name));
                try
                    fid = fopen(out_file_name, 'Wb');
                    
                    str_tmp = sprintf('# XYZ Position file generated on %s \n', now_time.toString('dd-mmm-yyyy HH:MM'));
                    str_tmp = sprintf('%s#--------------------------------------------------------------------------------\n', str_tmp);
                    str_tmp = sprintf('%s#LOCAL GEODETIC DATUM: WGS - 84\n\n', str_tmp);
                    str_tmp = sprintf('%s#  NUM | STATION NAME |     FIRST EPOCH |      LAST EPOCH |  N_POS\n', str_tmp);
                    for r = 1 : n_sta
                        str_tmp = sprintf('%s# %4d | %12s | %15s | %15s | %6d\n', str_tmp, r, sta_list(r).getMarkerName4Ch, datestr(st_time(r), 'yyyymmdd HHMMSS'), datestr(en_time(r), 'yyyymmdd HHMMSS'), out_list(r).getTimePositions.length());
                    end
                    str_tmp = sprintf('%s\n#--------------------------------------------------------------------------------\n', str_tmp);
                    str_tmp = sprintf('%s# STAR T OF POSITIONS \n', str_tmp);
                    str_tmp = sprintf('%s#--------------------------------------------------------------------------------\n\n', str_tmp);
                    if upper(type(1)) == 'X' % XYZ
                        str_tmp = sprintf('%s  NUM  STATION NAME      DATE   TIME          X (M)           Y (M)           Z (M)    FLAG\n\n', str_tmp);
                    elseif upper(type(1)) == 'E' % ENU
                        str_tmp = sprintf('%s  NUM  STATION NAME      DATE   TIME       East (M)       North (M)          Up (M)    FLAG\n\n', str_tmp);
                    else % if upper(type(1)) == 'G' % Geodetic
                        str_tmp = sprintf('%s  NUM  STATION NAME      DATE   TIME      lat (deg)       lon (deg)    Ortho. H (M)   Ellips. H (M)    FLAG\n\n', str_tmp);
                    end
                    
                    fprintf(fid, str_tmp);
                    
                    
                    n_rec = length(sta_list);
                    rf = Core.getReferenceFrame;
                    for r = 1 : n_rec
                        str_tmp = '';
                        if upper(type(1)) == 'X' % XYZ
                            coo = out_list(r).getPos;
                            coo = coo.getXYZ;
                        elseif upper(type(1)) == 'E' % ENU
                            coo = out_list(r).getPosENU;
                        else % if upper(type(1)) == 'G' % Geodetic
                            [lat, lon, h_ellips, h_ortho]= out_list(r).getPosGeodetic;
                        end
                        
                        t_pos = out_list(r).getTimePositions.getMatlabTime;
                        % Is the station fixed?
                        marker = sta_list(r).getMarkerName4Ch;
                        flag = iif(rf.getFlag(marker) == 2, 'F', 'P');
                        if upper(type(1)) == 'X' % XYZ
                            for e = 1 : size(coo, 1)
                                str_tmp = sprintf('%s %4d  %12s  %15s  %14.4f  %14.4f  %14.4f    %1s\n', str_tmp, r, marker, datestr(t_pos(e), 'yyyymmdd HHMMSS'), coo(e, 1), coo(e, 2), coo(e, 3), flag);
                            end
                        elseif upper(type(1)) == 'E' % ENU
                            for e = 1 : size(coo, 1)
                                str_tmp = sprintf('%s %4d  %12s  %15s  %14.4f  %14.4f  %14.4f    %1s\n', str_tmp, r, marker, datestr(t_pos(e), 'yyyymmdd HHMMSS'), coo(e, 1), coo(e, 2), coo(e, 3), flag);
                            end
                        else % if upper(type(1)) == 'E' % ENU
                            for e = 1 : size(lat, 1)
                                str_tmp = sprintf('%s %4d  %12s  %15s  %14.9f  %14.9f  %14.4f  %14.4f    %1s\n', str_tmp, r, marker, datestr(t_pos(e), 'yyyymmdd HHMMSS'), lat(e) / pi * 180, lon(e) / pi * 180, h_ortho(e), h_ellips(e), flag);
                            end
                        end
                        fprintf(fid, str_tmp);
                    end
                    fclose(fid);
                    Core.getLogger.addStatusOk(sprintf('Exporting completed successfully'));
                catch ex
                    Core.getLogger.addError(sprintf('Exporting failed'));
                    Core_Utils.printEx(ex);
                end
            end
        end

        function exportMat(sta_list)
            % Export the receiver into a MATLAB file (work properties is not saved )
            %
            % SYNTAX
            %   sta_list.exportMat()
            
            for r = 1 : numel(sta_list)
                try
                    % Get time span of the receiver
                    time = sta_list(r).getTime().getEpoch([1 sta_list(r).getTime().length()]);
                    
                    fname = fullfile(Core.getState.getOutDir(), sprintf('full_%s-%s-%s-rec%04d%s', sta_list(r).getMarkerName4Ch, time.first.toString('yyyymmdd_HHMMSS'), time.last.toString('yyyymmdd_HHMMSS'), r, '.mat'));
                    
                    rec = sta_list(r);
                    tmp_work = rec.work; % back-up current out
                    rec.work = Receiver_Work_Space(rec);
                    save(fname, 'rec');
                    rec.work = tmp_work;
                    
                    rec.log.addStatusOk(sprintf('Receiver %s: %s', rec.getMarkerName4Ch, fname));
                catch ex
                    sta_list(r).log.addError(sprintf('saving Receiver %s in matlab format failed: %s', sta_list(r).getMarkerName4Ch, ex.message));
                end
            end
        end
        
        function exportPlainMat(sta_list, out_file_name)
            % Export the Results of the processing as a .mat file
            % without using objects
            %
            % SYNTAX:
            %    sta_list.exportPlainMat(<out_file_name>)
            
            log = Core.getLogger;
            state = Core.getState;
            if nargin < 2 || isempty(out_file_name)
                out_dir = state.getOutDir();
                if numel(sta_list) ~= 1
                    % contains multiple stations
                    out_file_name = fullfile(out_dir, sprintf('plain_output_full_%s.mat',GPS_Time.now.toString('yyyymmdd_HHMMSS')));
                else
                    % contains just one station
                    out_file_name = fullfile(out_dir, sprintf('plain_output_%s_%s.mat', upper(sta_list(1).getMarkerName4Ch()), GPS_Time.now.toString('yyyymmdd_HHMMSS')));
                end
            else
                if sum(out_file_name == filesep()) == 0
                    out_dir = state.getOutDir();
                    out_file_name = fullfile(out_dir, out_file_name);
                end
            end
            log.addMarkedMessage(sprintf('Exporting output to "%s"',out_file_name));
            
            for r = 1 : numel(sta_list)                
                if ~sta_list(r).isEmptyOut_mr
                    rec(r).short_name_4ch = sta_list(r).getMarkerName4Ch();
                    log.addMessage(log.indent(sprintf(' - appending %s', rec(r).short_name_4ch)));
                    rec(r).description_short = sta_list(r).getMarkerName();
                    rec(r).a_priori_xyz = sta_list(r).getMedianPosXYZ();
                    if ~isempty(sta_list(r).work) && (~sta_list(r).work.isEmpty)
                        rec(r).active_constellations = unique(sta_list(r).work.system);
                    else
                        cc = core.getConstellationCollector;
                        rec(r).active_constellations = cc.getActiveSysChar();
                    end
                    rec(r).observer = sta_list(r).observer;
                    rec(r).agency = sta_list(r).agency;
                    rec(r).ant_type = sta_list(r).ant_type;
                    [n_sat, epoch_time] = sta_list(r).getNumSat();
                    rec(r).epoch_time = epoch_time.getUnixTime;
                    rec(r).n_sat = n_sat(:);
                    rec(r).ztd = sta_list(r).getZtd_mr();
                    rec(r).zwd = sta_list(r).getZwd_mr();
                    rec(r).pwv = sta_list(r).getPwv_mr();
                    [rec(r).pressure, rec(r).temperature, rec(r).humidity] = sta_list(r).getPTH_mr();
                    [pos_xyz, pos_time] = sta_list(r).getPosXYZ();
                    rec(r).pos_time = pos_time{1}.getUnixTime();
                    rec(r).pos_xyz = pos_xyz{1};
                    rec(r).rate = core.state.getSessionDuration();
                end
            end
            log.addMessage(log.indent('Writing the file...'));
            save(out_file_name, 'rec', '-v7');
            log.addStatusOk('Export completed successfully');
        end
        
        function exportSlantTD(sta_list, time_rate, mode, format, flag_work)
            % Export Slant Total Delay
            % 
            % INPUT 
            %   time_rate   it could be a rate in seconds - default 21600 (6 hours)
            %               or an object of class GPS_Time
            %               (e.g. GPS_Time.fromString(2013/06/04 12:00:00))
            %
            %   mode        indicate if it is necessary to add the residuals
            %               0 - do not add residuals
            %               1 - add residuals 
            %               2 - add spline smoothed residuals (1 spline every 900 seconds)
            %
            %   format      'compact' - MJD SAT SWD EL AZ (default)
            %               'extended' - extended - Year DOY SOD SWD SAT SAT_X SAT_Y SAT_Z Receiver REC_X REC_Y REC_Z elAngle(rad) elAngle(deg) Azimuth(rad) Azimuth(deg)
            %
            %   flag_work  use the receiver work_space instead of output
            %
            %
            % NOTE: at the moment adding residuals works only if the data have been porcessed 
            %       with the combined engine in iono-free mode
            %
            % SYNTAX
            %   this.exportSlantTD(this, time_rate, mode, format, flag_work)
            
            log = Core.getLogger();
            if nargin < 4 || isempty(format)
                format = 'compact';
            end
            if nargin < 3 || isempty(mode)
                mode = 0;
            end
            if nargin < 2 || isempty(time_rate)
                time_rate = 21600;
            end
            err_code = 0;
            if nargin >= 5 && ~isempty(flag_work) && flag_work
                % use receiver work_space
                sta_list = sta_list(~sta_list.isEmptyWork_mr);
                for s = 1 : numel(sta_list)
                    err_code = err_code + sta_list(s).work.exportSlantDelay(time_rate, mode, format, iif(format(1) == 'c' && s == 1, false, true), false);
                end
            else
                % use receiver out
                sta_list = sta_list(~sta_list.isEmptyOut_mr);
                for s = 1 : numel(sta_list)
                    err_code = err_code + sta_list(s).out.exportSlantDelay(time_rate, mode, format, iif(format(1) == 'c' && s == 1, false, true), false);
                end                
            end
            if err_code > 0
                log.addError('Export failed');
            else
                log.addStatusOk('Export completed');
            end
        end
        
        function exportSlantWD(sta_list, time_rate, mode, format, flag_work)
            % Export Slant Wet Delay
            % 
            % INPUT 
            %   time_rate   it could be a rate in seconds - default 21600 (6 hours)
            %               or an object of class GPS_Time
            %               (e.g. GPS_Time.fromString(2013/06/04 12:00:00))
            %
            %   mode        indicate if it is necessary to add the residuals
            %               0 - do not add residuals
            %               1 - add residuals 
            %               2 - add spline smoothed residuals (1 spline every 900 seconds)
            %
            %   format      'compact' - MJD SAT SWD EL AZ (default)
            %               'extended' - extended - Year DOY SOD SWD SAT SAT_X SAT_Y SAT_Z Receiver REC_X REC_Y REC_Z elAngle(rad) elAngle(deg) Azimuth(rad) Azimuth(deg)
            %
            %   flag_work  use the receiver work_space instead of output
            %
            %
            % NOTE: at the moment adding residuals works only if the data have been porcessed 
            %       with the combined engine in iono-free mode
            %
            % SYNTAX
            %   this.exportSlantWD(this, time_rate, mode, format, flag_work)
            
            log = Core.getLogger();
            if nargin < 4 || isempty(format)
                format = 'compact';
            end
            if nargin < 3 || isempty(mode)
                mode = 0;
            end
            if nargin < 2 || isempty(time_rate)
                time_rate = 21600;
            end
            err_code = 0;
            if nargin >= 5 && ~isempty(flag_work) && flag_work
                % use receiver work_space
                sta_list = sta_list(~sta_list.isEmptyWork_mr);
                for s = 1 : numel(sta_list)
                    err_code = err_code + sta_list(s).work.exportSlantDelay(time_rate, mode, format, iif(format(1) == 'c' && s == 1, false, true), true);
                end
            else
                % use receiver out
                sta_list = sta_list(~sta_list.isEmptyOut_mr);
                for s = 1 : numel(sta_list)
                    err_code = err_code + sta_list(s).out.exportSlantDelay(time_rate, mode, format, iif(format(1) == 'c' && s == 1, false, true), true);
                end                
            end
            if err_code > 0
                log.addError('Export failed');
            else
                log.addStatusOk('Export completed');
            end
        end
        
        function exportMultiPath(sta_list)
            % Export the current multi-path mitigation model stored in GNSS_Station
            % Each file is stored in a folder with name = marker 4ch
            %
            % SYNTAX:
            %    sta_list.exportMultiPath(<out_file_name>)
            
            log = Core.getLogger;
            state = Core.getState;
            flag_export = 0;
            out_dir = state.getMPDir();
            for r = 1 : numel(sta_list)
                try
                    rec = sta_list(r);
                    if ~isempty(rec.ant_mp)
                        if exist(fullfile(out_dir, rec.getMarkerName4Ch), 'dir') ~= 7
                            mkdir(fullfile(out_dir, rec.getMarkerName4Ch));
                        end
                        ant_mp = rec.ant_mp; %#ok<NASGU>
                        out_file_name = fullfile(out_dir, rec.getMarkerName4Ch, sprintf('%s_mp_%s_%s.mat', ...
                            rec.getMarkerName4Ch, ...
                            ant_mp.time_lim.first.toString('yyyymmddHHMM'), ...
                            ant_mp.time_lim.last.toString('yyyymmddHHMM')));
                        
                        log.addMarkedMessage(sprintf('Exporting multipath model to "%s"',out_file_name));
                        save(out_file_name, 'ant_mp');
                        flag_export = flag_export + 1;
                    end
                catch ex
                    Core_Utils.printEx(ex);
                end
            end
            if flag_export == 0 && numel(sta_list) > 0
                log.addError('Multipath mitigation model export failed, no models found!');
            else
                log.addStatusOk(sprintf('%d of %d multipath mitigation models exported successfully', flag_export, numel(sta_list)));
            end
        end
        
        function exportHydroNET(sta_list)
            % Export the troposphere into a hydroNet readable data format file
            % The data exported are:
            %  - ztd
            %  - zwd
            %  - east gradient
            %  - north gradient
            %  - time_utc in matlab format
            %
            % SYNTAX
            %   this.exportHydroNET
            try
                sta_list = sta_list(~sta_list.isEmptyOut_mr);
                if ~isempty(sta_list)
                    state = Core.getCurrentSettings;
                    min_time = state.sss_date_start.getCopy();
                    max_time = state.sss_date_stop.getCopy();
                    %                 for r = 1 : numel(sta_list)
                    %                     min_time.append(sta_list(r).out.time.minimum);
                    %                     max_time.append(sta_list(r).out.time.maximum);
                    %                 end
                    %                 min_time = min_time.minimum;
                    %                 max_time = max_time.maximum;
                    [year,doy] = min_time.getDOY();
                    t_start = min_time.toString('HHMM');
                    
                    out_dir = fullfile( Core.getState.getOutDir(), sprintf('%4d', year), sprintf('%03d',doy));
                    if ~exist(out_dir, 'file')
                        mkdir(out_dir);
                    end
                    
                    if length(sta_list) == 1
                        prefix = sta_list(1).getMarkerName4Ch;
                    else
                        prefix = 'HNTD';
                    end
                    fname = sprintf('%s',[out_dir filesep prefix sprintf('%04d%03d_%4s_%d', year, doy, t_start, round(max_time.last()-min_time.first())) 'HN.csv']);
                    fid = fopen(fname,'Wb');
                    
                    % Switch to UTC
                    min_time.toUtc();
                    max_time.toUtc();
                    
                    % write header
                    fprintf(fid,'[Variables]\n');
                    fprintf(fid,'Code,Name,Unit\n');
                    fprintf(fid,'ZTD,Zenith Total Delay,m\n');
                    fprintf(fid,'ZWD,Zenith Wet Delay,m\n');
                    fprintf(fid,'GE,East Gradient,m\n');
                    fprintf(fid,'GN,North Gradient,m\n');
                    fprintf(fid,'\n');
                    fprintf(fid,'[Locations]\n');
                    fprintf(fid,'Code,Name,X,Y,Z,EPSG\n'); % XYZ geocentric
                    %fprintf(fid,'Code, Name, latitude, longitude, h_ortho, EPSG\n');
                    for r = 1 : numel(sta_list)
                        coo = Coordinates.fromXYZ(sta_list(r).out.xyz(1,:));%Coordinates.fromXYZ(, sta_list(r).out.time_pos);
                        %[x,y,h_ellipse,zone] = coo.getENU();
                        %                     ondu = coo.getOrthometricCorrection();
                        %                     h_ortho = h_ellipse - ondu;
                        % if zone(4) < 'N'
                        %     hemi = '7';
                        % else
                        %     hemi = '6';
                        % end
                        
                        [lat, lon, ~, h_o] = coo.getGeodetic();
                        %fprintf(fid,'%s,%s,%0.4f,%0.4f,%0.4f,%s\n',sta_list(r).getMarkerName4Ch, sta_list(r).getMarkerName,x,y,z,'4326');  % XYZ geocentric
                        fprintf(fid,'%s,%s,%.8f,%.8f,%0.3f,%s\n', sta_list(r).getMarkerName4Ch, sta_list(r).getMarkerName, lon / pi * 180, lat / pi * 180, h_o, '4326');
                    end
                    fprintf(fid,'\n');
                    
                    fprintf(fid,'[Time]\n');
                    fprintf(fid,'UTC time zone offset,Model date,Start date,End date\n');
                    fprintf(fid,'|+0000,,%s,%s\n',min_time.toString('yyyy-mm-dd HH:MM:SS'),max_time.toString('yyyy-mm-dd HH:MM:SS'));
                    fprintf(fid,'\n');
                    
                    fprintf(fid,'[Data]\n');
                    fprintf(fid,'Date time,Location code,Variable code,Value,Quality,Availability\n');
                    for r = 1 : numel(sta_list)
                        if ~sta_list(r).isEmpty && ~sta_list(r).isEmptyOut_mr && ~isempty(sta_list(r).out.quality_info.s0) && max(sta_list(r).out.quality_info.s0) < 0.10
                            time = sta_list(r).out.getTime();
                            time.toUtc();
                            
                            ztd = sta_list(r).out.getZtd();
                            zwd = sta_list(r).out.getZwd();
                            [gn,ge ] =  sta_list(r).out.getGradient();
                            mk_code = sta_list(r).getMarkerName4Ch;
                            for t = 1 : time.length
                                time_str = time.getEpoch(t).toString('yyyy-mm-dd HH:MM:SS');
                                fprintf(fid,'%s,%s,ZTD,%0.4f,,1\n',time_str,mk_code,ztd(t));
                                fprintf(fid,'%s,%s,ZWD,%0.4f,,1\n',time_str,mk_code,zwd(t));
                                fprintf(fid,'%s,%s,GN,%0.5f,,1\n',time_str,mk_code,gn(t));
                                fprintf(fid,'%s,%s,GE,%0.5f,,1\n',time_str,mk_code,ge(t));
                            end
                        else
                            if isempty(max(sta_list(r).out.quality_info.s0))
                                Core.getLogger().addWarning(sprintf('no solution have been found, station skipped'));
                            else
                                Core.getLogger().addWarning(sprintf('s02 (%f m) too bad, station skipped', max(sta_out(r).out.quality_info.s0)));
                            end
                        end
                    end
                    fclose(fid);
                    log = Core.getLogger;
                    log.addStatusOk(sprintf('Tropo saved into: "%s"', fname));
                    log.addStatusOk('Export completed successfully');
                else
                    log = Core.getLogger;
                    log.addWarning(sprintf('no solution have been found, station skipped'));
                end
            catch ex
                if all(sta_list.isEmptyOut_mr)
                    log = Core.getLogger;
                    log.addWarning(sprintf('no solution have been found, station skipped'));
                else
                    Core_Utils.printEx(ex);
                    log = Core.getLogger;
                    log.addError(sprintf('saving Tropo in CSV (HydroNet) format failed: "%s"', ex.message));
                end
            end
        end
    end
    %% METHODS GETTER
    % ==================================================================================================================================================
    
    methods
        % standard utility
        function toString(sta_list)
            % Display on screen information about the receiver
            %
            % INPUT
            %   sta_list    list of receivers
            %
            % SYNTAX
            %   this.toString(sta_list);
            log = Core.getLogger();
            for i = 1:length(sta_list)
                if ~isempty(sta_list(i))
                    log.starSeparator();
                    log.addMarkedMessage(sprintf('Receiver %s\n Object created at %s', sta_list(i).getMarkerName(), sta_list(i).creation_time.toString));
                    log.starSeparator();

                    if ~sta_list(i).work.isEmpty()
                        sta_list(i).work.toString();
                    end
                    if ~sta_list(i).out.isEmpty()
                        sta_list(i).out.toString();
                    end
                end
            end
        end

        function ant_mp = getAntennaMultiPath(this)
            % Get the structure containing multipath information for the receiver
            %
            % SYNTAX 
            %   ant_mp = this.getAntennaMultiPath()
            
            % If no multipath map is loaded try to load the map
            log = Core.getLogger;
            if isempty(this.ant_mp)
                state = Core.getState;
                marker_name = upper(this.getMarkerName4Ch);
                out_dir = fullfile(state.getMPDir(), marker_name);
                out_file_name = fullfile(out_dir, sprintf('%s_mp_*.mat', marker_name));
                file_info = dir(out_file_name);
                if isempty(file_info)
                    marker_name = lower(this.getMarkerName4Ch);
                    out_dir = fullfile(state.getMPDir(), marker_name);
                    out_file_name = fullfile(out_dir, sprintf('%s_mp_*.mat', marker_name));
                    file_info = dir(out_file_name);
                elseif isempty(file_info)
                    marker_name = (this.getMarkerName4Ch);
                    out_dir = fullfile(state.getMPDir(), marker_name);
                    out_file_name = fullfile(out_dir, sprintf('%s_mp_*.mat', marker_name));
                    file_info = dir(out_file_name);
                end
                date = [];
                idf = [];
                for f = 1 : numel(file_info)
                    try
                        date_lim = regexp(file_info(f).name, '(?<=_)[0-9]*(?=_|.mat)', 'match');
                        date_lim = [datenum(date_lim{1}, 'yyyymmddHHMM'); datenum(date_lim{2}, 'yyyymmddHHMM')];
                        date = [date; date_lim'];
                        idf = [idf; f f];
                    catch ex
                        log.addWarning(sprintf('File not valid found in MP dir "%s"', file_info(f).name));
                    end
                end 
                % Find the closer multipath map file
                if ~isempty(date)
                    try
                        if isempty(this.work) || isempty(this.work.getTime)
                            sss_center = core.state.sss_date_start.getMatlabTime + core.state.sss_duration/86400 * (core.getCurrentSession -0.5);
                        else
                            sss_center = this.work.getTime.getCentralTime.getMatlabTime;
                        end
                        % check if the data is inside a MP map
                        closer_fid = idf(find((date(:,2) >= sss_center) & (date(:,1) <= sss_center), 1, 'last'), 2);
                        % Always take the newer map in case of overlap
                        if isempty(closer_fid)
                            % Check for the closer map
                            date = date(:);
                            idf = idf(:);
                            [date, id_sort] = sort(date);
                            [~, id_min] = min(abs(date - sss_center));
                            closer_fid = idf(id_sort(id_min));
                        end
                        file_name = fullfile(out_dir, file_info(closer_fid).name);
                        log.addMarkedMessage(sprintf('%s - Importing Multipath mitigation model from "%s"',  marker_name, file_name));
                        load(file_name, 'ant_mp');
                        this.ant_mp = ant_mp;
                    catch ex
                        log.addWarning(sprintf('Loading Multipath mitigation model for %s failed', marker_name));
                    end
                end
            end
            ant_mp = this.ant_mp;
        end
        
        function sys_c = getActiveSys(this)
            % Get the active system stored into the object
            %
            % SYNTAX 
            %   sys_c = this.getActiveSys()
            
            % Select only the systems still present in the object           
            if ~isempty(this.out) && ~(this.out.isEmpty) && ~isempty(this.out.used_sys_c)
                sys_c = this.out.used_sys_c;
            elseif ~isempty(this.work) && ~(this.work.isEmpty)
                sys_c = this.work.getActiveSys;
            else                
                sys_c = this.getCC.getActiveSysChar;
            end
        end
        
        function cc = getCC(this)
            % Get Constellation collector
            %
            % SYNTAX
            %   cc = this.getCC()
            cc = Core.getState.getConstellationCollector;
        end
        
        function id = getStationId(sta_list, marker_name)
            % Given a marker_name get the sequencial id of a station
            %
            % INPUT
            %   sta_list      list of receivers
            %   marker_name   4 letter marker name
            %
            % SYNTAX
            %   id = getStationId(this, marker_name)
            marker4ch_list = '';
            for r = 1 : numel(sta_list)
                try
                    marker4ch_list(r, :) = char(sta_list(r).getMarkerName4Ch);
                catch
                    % the name is shorter or missing => ignore
                end
            end
            id = find(Core_Utils.code4Char2Num(upper(marker4ch_list)) == Core_Utils.code4Char2Num(upper(marker_name)));
        end

        function [req_rec, id_rec] = get(sta_list, marker_name)
            % Get the receivers with a certain Marker name (case unsensitive)
            %
            % SYNTAX
            %   req_rec = sta_list.get(marker_name)
            %
            % SEE ALSO
            %   sta_list.printStationList()
            req_rec = [];
            id_rec = [];
            for r = 1 : size(sta_list,2)
                rec = sta_list(~sta_list(:,r).isEmpty_mr ,r);
                if not(rec.isEmpty)
                    if strcmpi(rec(1).getMarkerName, marker_name)
                        id_rec = [id_rec r]; %#ok<AGROW>
                        req_rec = [req_rec sta_list(:,r)]; %#ok<AGROW>
                    elseif strcmpi(rec(1).getMarkerName4Ch, marker_name)
                        id_rec = [id_rec r]; %#ok<AGROW>
                        req_rec = [req_rec sta_list(:,r)]; %#ok<AGROW>
                    end
                end
            end
        end

        function marker_name = getMarkerName(this)
            % Get the Marker name as specified in the RINEX file
            %
            % SYNTAX
            %   marker_name = getMarkerName(this)
            marker_name = this.marker_name;
            if isempty(marker_name)
                marker_name = this.getMarkerName4Ch();
            end
        end

        function printStationList(sta_list)
            % Print the list of station / markers
            %
            % SYNTAX
            %   sta_list.printStationList()
            log = Core.getLogger();
            if numel(sta_list) > 0
                log.addMessage('List of available stations:');
                for r = 1 : numel(sta_list)
                    try
                        log.addMessage(sprintf('%4d) %s - %s', r, char(sta_list(r).getMarkerName4Ch), sta_list(r).getMarkerName));
                    catch
                        % the name is shorter or missing => ignore
                    end
                end
            end
        end

        function marker_name = getMarkerName4Ch(sta_list)
            % Get the Marker name as specified in the file name
            % (first four characters)
            %
            % SYNTAX
            %   marker_name = getMarkerName4Ch(this)
            marker_name = cell(numel(sta_list), 1);
            for s = 1 : numel(sta_list)
                if ~isempty(sta_list(s).marker_name)
                    marker_name{s} = sta_list(s).marker_name;
                else
                    if isempty(sta_list(s).work.rinex_file_name)
                        marker_name{s} = 'NONE';
                    else
                        marker_name{s} = File_Name_Processor.getFileName(sta_list(s).work.rinex_file_name);
                    end
                end
                marker_name{s} = marker_name{s}(1 : min(4, length(marker_name{s})));
            end
            if numel(sta_list) == 1
                marker_name = marker_name{1};
            end
        end

        function out_prefix = getOutPrefix(this)
            % Get the name for exporting output (valid for dayly output)
            %   - marker name 4ch (from rinex file name)
            %   - 4 char year
            %   - 3 char doy
            %
            % SYNTAX
            %   out_prefix = this.getOutPrefix()
            if this.out.length == 0
                time = this.work.time.getCopy;
            else
                time = this.out.time.getCopy;
            end
            [year, doy] = time.getCentralTime.getDOY();
            out_prefix = sprintf('%s_%04d_%03d_', this.getMarkerName4Ch, year, doy);
        end

        function is_empty = isEmpty_mr(sta_list)
            % Return if the object does not cantain any observations (work) or results (out)
            %
            % SYNTAX
            %   is_empty = this.isEmpty_mr();
            %
            % SEE ALSO
            %   isEmptyOut_mr isEmptyWork_mr
            is_empty =  false(numel(sta_list), 1);
            for r = 1 : numel(sta_list)
                is_empty(r) =  sta_list(r).work.isEmpty() && sta_list(r).out.isEmpty();
            end
        end
        
        function has_phases = hasPhases_mr(sta_list)
            % Return if the object does not cantain any observations (work) or results (out)
            %
            % SYNTAX
            %   is_empty = this.isEmpty_mr();
            %
            % SEE ALSO
            %   isEmptyOut_mr isEmptyWork_mr
            has_phases = false(numel(sta_list), 1);
            for r = 1 : numel(sta_list)
                has_phases(r) =  ~sta_list(r).work.isEmpty() && sta_list(r).work.hasPhases();
            end
        end
        
         function has_phases = hasRangeObs_mr(sta_list)
            % Return if the object does not cantain any observations (work) or results (out)
            %
            % SYNTAX
            %   is_empty = this.hasRangeObs_mr();
            %
            % SEE ALSO
            %   isEmptyOut_mr isEmptyWork_mr
            has_phases = false(numel(sta_list), 1);
            for r = 1 : numel(sta_list)
                has_phases(r) =  ~sta_list(r).work.isEmpty() && sta_list(r).work.hasRangeObs();
            end
        end
        

        function is_empty = isEmptyWork_mr(sta_list)
            % Return if the object work does not cantains any observation
            %
            % SYNTAX
            %   is_empty = this.isEmptyWork_mr();
            %
            % SEE ALSO
            %   isEmpty_mr isEmptyOut_mr
            is_empty =  false(numel(sta_list), 1);
            for r = 1 : numel(sta_list)
                is_empty(r) =  sta_list(r).work.isEmpty();
            end
        end

        function is_empty = isEmptyOut_mr(sta_list)
            % Return if the object out does not cantains any results
            %
            % SYNTAX
            %   is_empty = this.isEmptyOut_mr();
            %
            % SEE ALSO
            %   isEmpty_mr isEmptyWork_mr
            is_empty =  false(numel(sta_list), 1);
            for r = 1 : numel(sta_list)
                is_empty(r) =  sta_list(r).out.isEmpty();
            end
        end

        function is_empty = isEmpty(this)
            % Return if the object does not cantains any observation
            %
            % SYNTAX
            %   is_empty = this.isEmpty();

            is_empty = isempty(this) || ((isempty(this.work) || this.work.isEmpty()) && (isempty(this.out) || this.out.isEmpty()));
        end

        function time = getTime(this)
            % return the time stored in the object out
            %
            % OUTPUT
            %   time     GPS_Time
            %
            % SYNTAX
            %   xyz = this.getTime()
            time = this.out.getTime();
        end

        function n_epo = getNumEpochs(sta_list)
            % Return the number of epochs stored in work
            %
            % SYNTAX
            %   len = this.getNumEpochs();
            n_epo =  zeros(numel(sta_list), 1);
            for r = 1 : numel(sta_list)
                n_epo(r) =  sta_list(r).work.time.length();
            end
            n_epo = sum(n_epo);
        end

        function n_sat = getMaxSat(sta_list, sys_c)
            % get the number of satellites stored in the object work
            %
            % SYNTAX
            %   n_sat = getNumSat(<sys_c>)
            n_sat = zeros(numel(sta_list),1);
            for r = 1 : size(sta_list, 2)
                rec(r) = sta_list(r);
                if nargin == 2
                    n_sat(r) = rec.work.getMaxSat(sys_c);
                elseif nargin == 1
                    n_sat(r) = rec.work.getMaxSat();
                end
            end
        end

        function n_sat = getMaxNumSat(sta_list, sys_c)
            % get the number of maximum theoretical satellites stored in the object
            %
            % SYNTAX
            %   n_sat = getMaxNumSat(<sys_c>)

            n_sat = zeros(numel(sta_list),1);
            for r = 1 : size(sta_list, 2)
                rec(r) = sta_list(r);
                if nargin == 2
                    n_sat(r) = rec.work.getMaxNumSat(sys_c);
                elseif nargin == 1
                    n_sat(r) = rec.work.getMaxNumSat();
                end
            end
        end

        function [time_lim_small, time_lim_large] = getWorkTimeSpan(this)
            % return a GPS_Time containing the first and last epoch stored in the Receiver
            %
            % OUTPUT
            %   time_lim_small     GPS_Time (first and last) epoch of the smaller interval
            %   time_lim_large     GPS_Time (first and last) epoch of the larger interval
            %
            % SYNTAX
            %   [time_lim_small, time_lim_large] = getWorkTimeSpan(this);
            %
            time_lim_small = this(1).work.time.first;
            tmp_small = this(1).work.time.last;
            time_lim_large = time_lim_small.getCopy;
            tmp_large = tmp_small.getCopy;
            for r = 2 : numel(this)
                if time_lim_small < this(r).work.time.first
                    time_lim_small = this(r).work.time.first;
                end
                if time_lim_large > this(r).work.time.first
                    time_lim_large = this(r).work.time.first;
                end

                if tmp_small > this(r).work.time.last
                    tmp_small = this(r).work.time.last;
                end
                if tmp_large < this(r).work.time.last
                    tmp_large = this(r).work.time.last;
                end
            end
            time_lim_small.append(tmp_small);
            time_lim_large.append(tmp_large);
        end

        function [time_lim_small, time_lim_large] = getOutTimeSpan(this)
            % return a GPS_Time containing the first and last epoch stored in the Receiver
            %
            % OUTPUT
            %   time_lim_small     GPS_Time (first and last) epoch of the smaller interval
            %   time_lim_large     GPS_Time (first and last) epoch of the larger interval
            %
            % SYNTAX
            %   [time_lim_small, time_lim_large] = getOutTimeSpan(this);
            %
            time_lim_small = this(1).out.time.first;
            tmp_small = this(1).out.time.last;
            time_lim_large = time_lim_small.getCopy;
            tmp_large = tmp_small.getCopy;
            for r = 2 : numel(this)
                if time_lim_small < this(r).out.time.first
                    time_lim_small = this(r).out.time.first;
                end
                if time_lim_large > this(r).out.time.first
                    time_lim_large = this(r).out.time.first;
                end

                if tmp_small > this(r).out.time.last
                    tmp_small = this(r).out.time.last;
                end
                if tmp_large < this(r).out.time.last
                    tmp_large = this(r).out.time.last;
                end
            end
            time_lim_small.append(tmp_small);
            time_lim_large.append(tmp_large);
        end

        function [rate] = getRate(this)
            % Get the rate of the output (or work if out is empty)
            %
            % SYNTAX
            %   rate = this.getRate();
            try
                if ~(isempty(this.out) || this.out.isEmpty)
                    rate = this.out.getTime.getRate;
                else
                    rate = nan;
                end
                if isnan(rate)
                    rate = this.work.getTime.getRate;
                end
            catch
                % if anything happen probably
                rate = nan;
            end
        end

        function [coo] = getPos(sta_list, flag_add_coo)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   coo     object array [ Coordinates ]
            %   
            % SYNTAX
            %   [coo] = sta_list.getPos(flag_add_coo)
            
            coo = Coordinates;
            if nargin == 1 || isempty(flag_add_coo) || flag_add_coo == 0
                for r = 1 : numel(sta_list)
                    coo(r) = sta_list(r).out.getPos();
                    if coo(r).time.isEmpty
                        coo(r).setTime(sta_list(r).out.getPositionTime());
                    end
                end
            else
                for r = 1 : numel(sta_list)
                    try
                        coo(r) = sta_list(r).out.add_coo(flag_add_coo).coo;
                    catch ex
                        % no coordinate found, skip
                    end
                end
            end
        end
        
        function printCrd_goGet(sta_list)
            % return at screen the CRD in goGet format
            %
            % SYNTAX
            %   coo = sta_list.printCrd()            
            crd = reshape(sta_list.getMedianPosXYZ, 3, numel(sta_list))';
            coo = sta_list.getPos;
            for r = 1 : numel(sta_list)
                [lat, lon] = coo(r).getGeodetic();
                lat = median(lat, 'omitnan');
                lon = median(lon, 'omitnan');
                fprintf('%4s_%04d_%04d     %13.4f %13.4f %13.4f 1\n', upper(sta_list(r).getMarkerName4Ch), round((90 - lat/pi*180)*10), round(mod(lon/pi*180,360)*10), crd(r, 1), crd(r, 2), crd(r, 3))
            end
        end

        function [xyz, time] = getPosXYZ(sta_list)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   xyz     XYZ coordinates cell array
            %
            % SYNTAX
            %   xyz = this.getPosENU()
            xyz = {};
            for r = 1 : numel(sta_list)
                xyz{r} = sta_list(r).out.getPosXYZ();
                time{r} = sta_list(r).out.getTimePositions();
            end
        end

        function [xyz, p_time, sta_ok] = getPosXYZ_mr(sta_list)
            % return the positions computed for the receiver
            % multi_rec mode (synced)
            %
            % OUTPUT
            %   xyz     XYZ coordinates synced matrix (n_epoch, 3, n_rec)
            %
            % SYNTAX
            %   [xyz, p_time, sta_ok] = sta_list.getPosXYZ_mr()

            sta_ok = find(~sta_list.isEmptyOut_mr());
            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list(sta_ok), [], true);

            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);
            
            n_rec = numel(sta_list);
            xyz = nan(size(id_sync, 1), 3, n_rec);
            for r = 1 : numel(sta_ok)                
                xyz_rec = sta_list(sta_ok(r)).out.getPosXYZ();
                id_rec = id_sync(:,r);
                xyz(~isnan(id_rec), :, sta_ok(r)) = xyz_rec(id_rec(~isnan(id_rec)), :);
            end
        end

        function [dist_3d, xyz_dist] = getDistFrom(sta_list, rec_ref)
            % Get distance [m] from the reference station rec_ref
            %
            % SYNTAX:
            %   [dist_3d, xyz_dist] = getDistFrom(this, rec_ref)
            
            xyz = zero2nan(sta_list.getMedianPosXYZ);
            xyz_dist = bsxfun(@minus, xyz, rec_ref.getMedianPosXYZ);
            dist_3d = sqrt(sum(xyz_dist.^2, 2));
        end

        function enu = getPosENU(sta_list)
            % return the positions computed for the receiver
            %
            % OUTPUT
            %   enu     enu coordinates cell array
            %
            % SYNTAX
            %   enu = sta_list.getPosENU()
            for r = 1 : numel(sta_list)
                enu{r} = sta_list(r).out.getPosENU();
            end
        end

        function [enu, p_time, sta_ok] = getPosENU_mr(sta_list)
            % return the positions computed for n receivers
            % multi_rec mode (synced)
            %
            % OUTPUT
            %   enu     enu synced coordinates
            %
            % SYNTAX
            %   [enu, p_time, sta_ok] = sta_list.getPosENU_mr()
            
            sta_ok = find(~sta_list.isEmptyOut_mr());
            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list(sta_ok), [], true);

            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);
            
            n_rec = numel(sta_list);
            enu = nan(size(id_sync, 1), 3, n_rec);
            for r = 1 : numel(sta_ok)                
                enu_rec = sta_list(sta_ok(r)).out.getPosENU();
                id_rec = id_sync(:,r);
                enu(~isnan(id_rec), :, sta_ok(r)) = enu_rec(id_rec(~isnan(id_rec)), :);
            end            
        end

        function xyz = getMedianPosXYZ(this)
            % return the computed median position of the receiver
            % MultiRec: works on an array of receivers
            %
            % OUTPUT
            %   xyz     geocentric coordinates
            %
            % SYNTAX
            %   xyz = this.getMedianPosXYZ()

            xyz = [];
            for r = 1 : numel(this)
                if isempty(median(this(r).out.getPosXYZ(), 1))
                    xyz = [xyz; nan(1,3)]; %#ok<AGROW>
                else
                    xyz = [xyz; median(this(r).out.getPosXYZ(), 1)]; %#ok<AGROW>
                end
            end
        end
        
        function enu = getMedianPosENU(sta_list)
            % return the computed median position of the receiver
            % MultiRec: works on an array of receivers
            %
            % OUTPUT
            %   enu         UTM East North Up coordinates
            %
            % SYNTAX
            %   enu = sta_list.getMedianPosENU();
            
            enu = nan(numel(sta_list), 3);
            for r = 1 : numel(sta_list)
                if sta_list(1).static
                    enu(r,:) = sta_list(r).out.getMedianPosENU();
                else
                    enu(r,:) = sta_list(r).out.getMedianPosENU();
                end
            end
        end
        
        function [lat, lon, h_ellips, h_ortho] = getMedianPosGeodetic(sta_list)
            % return the computed median position of the receiver
            % MultiRec: works on an array of receivers
            %
            % OUTPUT
            %   lat         latitude  [deg]
            %   lon         longitude [deg]
            %   h_ellips    ellipsoidical heigth [m]
            %   h_ortho     orthometric heigth [m]
            %
            % SYNTAX
            %   [lat, lon, h_ellips, h_ortho] = sta_list.getMedianPosGeodetic();

            lat = nan(numel(sta_list), 1);
            lon = nan(numel(sta_list), 1);
            h_ellips = nan(numel(sta_list), 1);
            h_ortho = nan(numel(sta_list), 1);
            for r = 1 : numel(sta_list)                
                if ~sta_list(r).isEmptyOut_mr
                    [lat_tmp, lon_tmp, h_ellips_tmp, h_ortho_tmp] = sta_list(r).out.getMedianPosGeodetic;
                else
                    [lat_tmp, lon_tmp, h_ellips_tmp, h_ortho_tmp] = sta_list(r).work.getMedianPosGeodetic;
                end
                if isnan(lat_tmp)
                    [lat_tmp, lon_tmp, h_ellips_tmp, h_ortho_tmp] = sta_list(r).work.getMedianPosGeodetic;
                end
                if sta_list(1).static
                    lat(r) = lat_tmp;
                    lon(r) = lon_tmp;
                    h_ellips(r) = h_ellips_tmp;
                    h_ortho(r) = h_ortho_tmp;
                else
                    lat{r} = lat_tmp;
                    lon{r} = lon_tmp;
                    h_ellips{r} = h_ellips_tmp;
                    h_ortho{r} = h_ortho_tmp;
                end
                
            end
        end

        function getChalmersString(sta_list)
            % Get the string of the station to be used in http://holt.oso.chalmers.se/loading/
            %
            % SYNTAX
            %   this.getChalmersString(sta_list);

            sta_list(1).log.addMarkedMessage('Chalmers ocean loading computation must be required manually:');
            sta_list(1).log.addMessage(sta_list(1).log.indent('go to http://holt.oso.chalmers.se/loading/ and request a BLQ file'));
            sta_list(1).log.addMessage(sta_list(1).log.indent('using ocean tide model FES2004'));
            % sta_list(1).log.addMessage(sta_list(1).log.indent('select also to compensate the values for the motion'));
            sta_list(1).log.addMessage(sta_list(1).log.indent('Use the following string for the station locations:'));
            sta_list(1).log.addMessage([char(8) '//------------------------------------------------------------------------']);

            for r = 1 : size(sta_list, 2)
                rec = sta_list(~sta_list(:,r).isEmpty, r);
                if ~isempty(rec)
                    xyz = rec.out.getMedianPosXYZ();
                    if isempty(xyz)
                        xyz = rec.work.getMedianPosXYZ();
                    end
                    sta_list(1).log.addMessage([char(8) sprintf('%-24s %16.4f%16.4f%16.4f', rec(1).getMarkerName4Ch, xyz(1), xyz(2),xyz(3))]);
                end
            end

            sta_list(1).log.addMessage([char(8) '//------------------------------------------------------------------------']);
        end

        function [pressure, temperature, humidity, p_time, id_sync] = getPTH_mr(sta_list)
            % Get synced data of TPH
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [pressure, temperaure, humidiy, p_time, id_sync] = this.getPTH_mr()

            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list);

            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);

            n_rec = numel(sta_list);
            pressure = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                id_rec(id_rec > length(sta_list(r).out.pressure)) = nan;
                pressure(~isnan(id_rec), r) = sta_list(r).out.pressure(id_rec(~isnan(id_rec)));
            end

            n_rec = numel(sta_list);
            temperature = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                id_rec(id_rec > length(sta_list(r).out.temperature)) = nan;
                temperature(~isnan(id_rec), r) = sta_list(r).out.temperature(id_rec(~isnan(id_rec)));
            end

            n_rec = numel(sta_list);
            humidity = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                id_rec(id_rec > length(sta_list(r).out.humidity)) = nan;
                humidity(~isnan(id_rec), r) = sta_list(r).out.humidity(id_rec(~isnan(id_rec)));
            end
        end
        
        
        function ztd_correction = getZtdReduction(sta_list, degree, xyh, flag_spatial)
            % Get the reduction due to height for ZTD
            % Estimate it with a quadratic 2d function on xy + 1d quadratic on height
            %
            % f(x,y,z) = a0 + a1 z + a2 z^2 + b1 x + b2 x^2 + b3 y + b4 y^2 + b5 * xy
            %
            % INPUT
            %   degree  height quadratic polynomial value
            %   xyh     [lon * cos(lat), lat, h_ortho] (angles in radians)
            %   flag spatial enable b coefficients
            %
            % SYNTAX
            %   ztd_correction = sta_list.getZtdReduction(degree, xyh, flag_spatial)
            
            if numel(sta_list) < 3
                ztd_correction = zeros(numel(sta_list), 1);
                log = Core.getLogger();
                log.addWarning('I cannot estimate an height correction with less than 3 stations');
                
                % Add ZHD here
            else
                data_med = median(sta_list.getZtd_mr, 'omitnan')';
                coo = Coordinates.fromXYZ(sta_list.getMedianPosXYZ());
                [lat, lon, ~, h_o] = coo.getGeodetic;
                if nargin < 4 || isempty(flag_spatial)
                    flag_spatial = false;
                end
                if flag_spatial
                    xyz = [lon .* cos(lat), lat, h_o];
                    if nargin < 3 || isempty(xyh)
                        xyh = xyz;
                    end
                    ztd_correction = Core_Utils.interp22nLS(xyz, data_med, degree, xyh);
                else
                    if nargin < 3 || isempty(xyh)
                        xyh = h_o;
                    else
                        xyh = xyh(:,3);
                    end
                    ztd_correction = Core_Utils.interp1LS(h_o, data_med, degree, xyh);
                end
                %try
                %    figure; subplot(2,1,1); plot(h_o, data_med, '.', lat_lon_h(:, 3), ztd_correction, '.'); subplot(2,1,2); plot(h_o, data_med - ztd_correction, '.');
                %catch
                %end
            end
        end
        
        
        function [ztd_res, p_time, ztd_height] = getReducedZtd_mr(sta_list, degree)
            % Reduce the ZTD of all the stations removing the component dependent with the altitude
            % Return synced ZTD
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [ztd_res, p_time, ztd_height] = sta_list.getReducedZtd_mr()

            med_ztd = median(sta_list.getZtd_mr, 'omitnan')';
            if nargin == 1
                degree = 2;
            end
            [~, ~, ~, h_o] = Coordinates.fromXYZ(sta_list.getMedianPosXYZ()).getGeodetic;

            ztd_height = Core_Utils.interp1LS(h_o, med_ztd, degree);
            [ztd, p_time] = sta_list.getZtd_mr();
            ztd_res = bsxfun(@minus, ztd', ztd_height)';
        end

        function [zwd_res, p_time, zwd_height] = getReducedZwd_mr(sta_list, degree)
            % Reduce the ZWD of all the stations removing the component dependent with the altitude
            % Return synced ZWD
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [ztd_res, p_time, ztd_height] = sta_list.getReducedZwd_mr()

            med_zwd = median(sta_list.getZwd_mr, 'omitnan')';
            if nargin == 1
                degree = 2;
            end
            [~, ~, ~, h_o] = Coordinates.fromXYZ(sta_list.getMedianPosXYZ()).getGeodetic;

            zwd_height = Core_Utils.interp1LS(h_o, med_zwd, degree);
            [zwd, p_time] = sta_list.getZwd_mr();
            zwd_res = bsxfun(@minus, zwd', zwd_height)';
        end
        
        function [pwv, p_time, id_sync] = getPwv_mr(sta_list)
            % Get synced data of pwv
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [pwv, p_time, id_sync] = sta_list.getPwv_mr()
            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list);

            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);

            n_rec = numel(sta_list);
            pwv = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                id_rec(id_rec > length(sta_list(r).out.pwv)) = nan;
                pwv(~isnan(id_rec), r) = sta_list(r).out.pwv(id_rec(~isnan(id_rec)));
            end
        end

        function [di, p_time, id_sync] = getSlantDispersionIndex_mr(sta_list)
            % Get synced data of dispersion index
            % Dispersion Index is computed as moving mean of var of slant total delay
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [zwd, p_time, id_sync] = this.getZwd_mr()
            %  [zwd, p_time, id_sync, tge, tgn] = this.getZwd_mr()
            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list);

            % Supposing all the station with the same constellation
            di = zeros(size(id_sync, 1), numel(sta_list));
            for r = 1 : numel(sta_list)
                try
                    sztd = sta_list(r).out.getSlantZTD(sta_list(r).slant_filter_win);
                    sztd = bsxfun(@minus, sztd, sta_list(r).out.ztd(id_sync(~isnan(id_sync(:, r)), r)));
                    di_tmp = 1e3 * sqrt(movmean(var(sztd', 'omitnan'), 1800 / sta_list(r).out.time.getRate, 'omitnan'));
                    di(~isnan(id_sync(:, r)), r) = di_tmp(id_sync(~isnan(id_sync(:, r)), r));
                catch
                    % missing station or invalid dataset
                end
            end
        end

        function [ztd_res, p_time, id_sync, ztd_height] = getZtdRes_mr(sta_list, degree)
            % Get synced data of ztd reduced by height effect
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [ztd, p_time, id_sync, ztd_height] = this.getZtd_mr()

            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list);

            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);

            n_rec = numel(sta_list);
            ztd_res = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                id_rec(id_rec > length(sta_list(r).out.ztd)) = nan;
                ztd_res(~isnan(id_rec), r) = sta_list(r).out.ztd(id_rec(~isnan(id_rec)));
            end
            
            med_ztd = median(ztd_res, 'omitnan')';
            [~, ~, ~, h_o] = Coordinates.fromXYZ(sta_list.getMedianPosXYZ()).getGeodetic;
            if nargin < 2 || ~isempty(degree)
                degree = 2;
            end
            ztd_height = Core_Utils.interp1LS(h_o, med_ztd, degree);
            [ztd_res, p_time] = sta_list.getZtd_mr();
            ztd_res = bsxfun(@minus, ztd_res', ztd_height)';
        end
        
        function [ztd, p_time, id_sync, tge, tgn] = getZtd_mr(sta_list)
            % Get synced data of ztd
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [ztd, p_time, id_sync] = this.getZtd_mr()
            %  [ztd, p_time, id_sync, tge, tgn] = this.getZtd_mr()

            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list);

            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);

            n_rec = numel(sta_list);
            ztd = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                id_rec(id_rec > length(sta_list(r).out.ztd)) = nan;
                ztd(~isnan(id_rec), r) = sta_list(r).out.ztd(id_rec(~isnan(id_rec)));
            end

            if nargout > 4
                tge = nan(size(id_sync));
                tgn = nan(size(id_sync));
                for r = 1 : n_rec
                    id_rec = id_sync(:,r);
                    id_rec(id_rec > length(sta_list(r).out.ztd)) = nan;
                    tge(~isnan(id_rec), r) = sta_list(r).out.tge(id_rec(~isnan(id_rec)));
                    tgn(~isnan(id_rec), r) = sta_list(r).out.tgn(id_rec(~isnan(id_rec)));
                end
            end
        end

        function [zwd_res, p_time, id_sync, zwd_height] = getZwdRes_mr(sta_list, degree)
            % Get synced data of zwd reduced by height effect
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [zwd, p_time, id_sync] = this.getZwd_mr()
            %  [zwd, p_time, id_sync, tge, tgn] = this.getZwd_mr()

            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list);

            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);

            n_rec = numel(sta_list);
            zwd_res = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                id_rec(id_rec > length(sta_list(r).out.zwd)) = nan;
                zwd_res(~isnan(id_rec), r) = sta_list(r).out.zwd(id_rec(~isnan(id_rec)));
            end
            
            med_zwd = median(zwd_res, 'omitnan')';
            [~, ~, ~, h_o] = Coordinates.fromXYZ(sta_list.getMedianPosXYZ()).getGeodetic;
            if nargin < 2 || ~isempty(degree)
                degree = 2;
            end
            zwd_height = Core_Utils.interp1LS(h_o, med_zwd, degree);
            [zwd_res, p_time] = sta_list.getZwd_mr();
            zwd_res = bsxfun(@minus, zwd_res', zwd_height)';
        end
            
        function [zwd, p_time, id_sync, tge, tgn] = getZwd_mr(sta_list)
            % Get synced data of zwd
            % MultiRec: works on an array of receivers
            %
            % SYNTAX
            %  [zwd, p_time, id_sync] = this.getZwd_mr()
            %  [zwd, p_time, id_sync, tge, tgn] = this.getZwd_mr()

            [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(sta_list);

            id_ok = any(~isnan(id_sync),2);
            id_sync = id_sync(id_ok, :);
            p_time = p_time.getEpoch(id_ok);

            n_rec = numel(sta_list);
            zwd = nan(size(id_sync));
            for r = 1 : n_rec
                id_rec = id_sync(:,r);
                id_rec(id_rec > length(sta_list(r).out.zwd)) = nan;
                zwd(~isnan(id_rec), r) = sta_list(r).out.zwd(id_rec(~isnan(id_rec)));
            end

            if nargout == 5
                tge = nan(size(id_sync));
                tgn = nan(size(id_sync));
                for r = 1 : n_rec
                    id_rec = id_sync(:,r);
                    id_rec(id_rec > length(sta_list(r).out.zwd)) = nan;
                    tge(~isnan(id_rec), r) = sta_list(r).out.tge(id_rec(~isnan(id_rec)));
                    tgn(~isnan(id_rec), r) = sta_list(r).out.tgn(id_rec(~isnan(id_rec)));
                end
            end
        end
        
        function [tropo, time, id_ko_cell] = getTropoPar(sta_list, par_name)
            % Get a tropo parameter among 'ztd', 'zwd', 'pwv', 'zhd'
            % Generic function multi parameter getter
            %
            % SYNTAX
            %  [tropo, p_time, id_ko] = sta_list.getAprZhd(sta_list, par_name)

            tropo = {};
            time = {};
            t_ref = 0;
            full_time = [];
            for r = 1 : numel(sta_list)
                tmp = sta_list(r).out.getTime();
                time{r} = tmp;
                if tmp.isEmpty
                    tropo{r} = [];
                else
                    if t_ref == 0
                        t_ref = round(time{r}.getCentralTime.getMatlabTime);
                    end
                    full_time = unique([full_time; round(time{r}.getRefTime(t_ref), 5)]);
                    
                    switch lower(par_name)
                        case 'ztd'
                            [tropo{r}] = sta_list(r).out.getZtd();
                        case 'zwd'
                            [tropo{r}] = sta_list(r).out.getZwd();
                            if isempty(tropo{r}) || all(isnan(zero2nan(tropo{r})))
                                [tropo{r}] = sta_list(r).out.getAprZwd();
                            end
                        case 'gn'
                            [tropo{r}] = sta_list(r).out.getGradient();
                        case 'ge'
                            [~,tropo{r}] = sta_list(r).out.getGradient();
                        case 'pwv'
                            [tropo{r}] = sta_list(r).out.getPwv();
                        case 'zhd'
                            [tropo{r}] = sta_list(r).out.getAprZhd();
                            if ~any(~isnan(tropo{r}))
                                ztd = sta_list(r).out.getZtd();
                                zwd = sta_list(r).out.getZwd();
                                tropo{r} = ztd-zwd;
                            end                                
                        case 'nsat'
                            [tropo{r}] = sta_list(r).out.getNSat();
                    end
                end
            end
            
            if numel(tropo) == 1
                tropo = tropo{1};
                time = time{1};
            end
            
            id_ko = [];
            id_ko_cell = {};
            if nargout == 3 && ~strcmpi(par_name, 'nsat') && ~strcmpi(par_name, 'zhd') && numel(sta_list) > 2
                Core.getLogger.addMessage('Compute outlier detection');
                tropo_sync = nan(numel(full_time), numel(sta_list));
                for r = 1 : numel(sta_list)
                    [~, id] = intersect(full_time, round(time{r}.getRefTime(t_ref), 5));
                    tropo_sync(id, r) = tropo{r}*1e2;
                end
                tropo_sync = bsxfun(@minus, tropo_sync, median(tropo_sync, 'omitnan'));
                tropo_sync = bsxfun(@minus, tropo_sync, median(tropo_sync, 2, 'omitnan'));
                id_ko = Core_Utils.snoopGatt(tropo_sync, 10, 4);
                for r = 1 : numel(sta_list)
                    [~, id] = intersect(full_time, round(time{r}.getRefTime(t_ref), 5));
                    id_ko_cell{r} = id_ko(id, r);
                end                
            end
        end

        function [tropo, time] = getZtd(sta_list)
            % Get ZTD
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getZtd()
            [tropo, time] = sta_list.getTropoPar('ztd');
        end

        function [n_sat, time] = getNumSat(sta_list)
            % Get the number of satellite in view per epoch
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getNumSat()
            [n_sat, time] = sta_list.getTropoPar('nsat');
        end

        function [tropo, time] = getZwd(sta_list)
            % Get ZWD
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getZwd()
            [tropo, time] = sta_list.getTropoPar('zwd');
        end

        function [tropo, time] = getPwv(sta_list)
            % Get PWV
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getPwv()
            [tropo, time] = sta_list.getTropoPar('pwv');
        end

        function [tropo, time] = getAprZhd(sta_list)
            % Get a-priori ZHD
            %
            % SYNTAX
            %  [tropo, p_time] = sta_list.getAprZhd()
            [tropo, time] = sta_list.getTropoPar('zhd');
        end
                        
        function [tropo_grid, x_grid, y_grid, time, tropo_height_correction, tropo_clim] = getTropoMap(sta_list, par_name, rate, flag_show)
            % Get interpolated map of tropospheric parameter
            % Resolution is determined by the dtm in use (2 * 0.029 degrees)
            % The map is computer only on ground (> 10m) + 1 degree of margin
            %
            % INPUT
            %   par_name    type of tropospheric parameter:
            %                - ztd
            %                - zwd
            %                - pwv
            %   rate        rate in seconds, nearest to closest observation 
            %               it should be a subsample of the data rate (e.g. 300 with 30s data)
            %   flag_show   if true show debug images
            %
            % SYNTAX
            %   [tropo_grid, x_grid, y_grid, time, tropo_height_correction] = sta_list.getTropoMap(par_name, rate)
            
            % Defining interpolation
            method = 'natural';
            dtm_size = 600; % keep dtm points smaller than dtm_size x dtm_size
            fun = @(dist) 0.2 * exp(-(dist)*1e1) + 0*exp(-(dist*5e1).^2);
            %fun = @(dist) 1./(dist+1e-5);

            if nargin < 4
                flag_show = false;
            end
            
            sta_list = sta_list(~sta_list.isEmptyOut_mr);
            switch lower(par_name)
                case 'ztd'
                    [tropo, s_time] = sta_list.getZtd_mr();
                case 'zwd'
                    [tropo, s_time] = sta_list.getZwd_mr();
                case 'dztd'
                    [tropo, s_time] = sta_list.getZtd_mr();
                    s_time = s_time.getEpoch(2:s_time.length);
                    s_time = s_time.addIntSeconds(s_time.getRate / 2);
                    tropo = (diff(tropo, 1)*1e2) / s_time.getRate * 3600;
                case 'dzwd'
                    [tropo, s_time] = sta_list.getZwd_mr();
                    s_time = s_time.getEpoch(2:s_time.length);
                    s_time = s_time.addIntSeconds(s_time.getRate / 2);
                    tropo = (diff(tropo, 1)*1e2) / s_time.getRate * 3600;
                case 'gn'
                    [~, s_time, ~, ~, tropo] = sta_list.getZwd_mr();
                case 'ge'
                    [~, s_time, ~, tropo, ~]  = sta_list.getZwd_mr();
                case 'dir'
                    [~, s_time, ~, ge, gn]  = sta_list.getZwd_mr();
                    tropo = atan2d(gn, ge) - 90;
                case 'pwv'
                    [tropo, s_time] = sta_list.getPwv_mr();
                case 'zhd'
                case 'nsat'
            end
            tropo = tropo * 1e2;
            
            if (nargin < 3) || isempty(rate)
                rate = 300; % 5 min
            end
            if numel(rate) == 1
                ss_rate = round(rate / s_time.getRate);
                time = s_time.getEpoch(1 : ss_rate : s_time.length());
            else
                time = s_time.getEpoch(rate);
            end

            med_tropo = mean(tropo, 1, 'omitnan')';

            degree = 2;
            coo = Coordinates.fromXYZ(sta_list.getMedianPosXYZ);
            [lat, lon, up, h_o] = coo.getGeodetic;
                        
            xyu = [lon / pi * 180, lat / pi * 180, up];
            
            % Remove missing stations
            id_ok = ~isnan(tropo'); % logical matrix of valid tropo
            
            % Generate interpolation grid
            x_span = max(xyu(:,1)) - min(xyu(:,1));
            y_span = max(xyu(:,2)) - min(xyu(:,2));
            border = 5;
            step = 0.05; % Grid step in degrees (only for defining limits)
            x_lim = [((floor(min(xyu(:,1)) / step) - border) * step), ((ceil(max(xyu(:,1)) / step) + border) * step)];
            y_lim = [((floor(min(xyu(:,2)) / step) - border) * step), ((ceil(max(xyu(:,2)) / step) + border) * step)];
            
            % Retrieve DTM model
            nwse = [max(y_lim), min(x_lim), min(y_lim), max(x_lim)];
            [dtm, lat, lon] = Core.getRefDTM(nwse, 'ortho', 'low');
            
            % I don't really need full resolution (maybe in a future)
            k = ceil(sqrt(numel(dtm)) / dtm_size);
            %k = 1;
            dtm = dtm(1:k:end, 1:k:end);
            lat = lat(1:k:end);
            lon = lon(1:k:end);
            
            % Correct data for height
            tropo_height = Core_Utils.interp1LS(h_o, med_tropo, degree);
            tropo_res = bsxfun(@minus, tropo', tropo_height)';            

            % Redifine the grid of the map to produce on the basis of the available DTM
            x_grid = lon';
            y_grid = flipud(lat)';
            x_step = median(diff(x_grid));
            y_step = median(diff(y_grid));
            [x_mg, y_mg] = meshgrid(x_grid, y_grid);

            % Generate computation mask (to limit points to interpolate)
            x_id = round((xyu(:,1) - x_grid(1)) / x_step) + 1;
            y_id = round((xyu(:,2) - y_grid(1)) / y_step) + 1;
            mask = zeros(size(y_grid, 2), size(x_grid, 2));
            mask((x_id - 1) * size(y_grid, 2) + y_id) = 1;
            % Refine mask keeping only points above sea level close to the station to process
            [xg, yg] = meshgrid(x_id, y_id); 
            % for Japan max was 25 -> elsewhere the bases are further away
            d = min(100, round(perc(noNaN(zero2nan(lower(sqrt((xg - xg').^2 + (yg - yg').^2)))), 0.15))); % 35% of min distance is used to enlarge the area of interpolation
            mask = (circConv2(mask, d) > 0) & (dtm >= -10);
            conv_mask = [0 0 1 1 1 0 0; ...
                         0 1 1 1 1 1 0; ...
                         1 1 1 1 1 1 1; ...
                         1 1 1 1 1 1 1; ...
                         1 1 1 1 1 1 1; ...
                         0 1 1 1 1 1 0; ...
                         0 0 1 1 1 0 0];
            conv_mask2 = [0 0 0 1 1 1 1 1 0 0 0; ...
                         0 0 1 1 1 1 1 1 1 0 0; ...
                         0 1 1 1 1 1 1 1 1 1 0; ...
                         1 1 1 1 1 1 1 1 1 1 1; ...
                         1 1 1 1 1 1 1 1 1 1 1; ...
                         1 1 1 1 1 1 1 1 1 1 1; ...
                         1 1 1 1 1 1 1 1 1 1 1; ...
                         1 1 1 1 1 1 1 1 1 1 1; ...
                         0 1 1 1 1 1 1 1 1 1 0; ...
                         0 0 1 1 1 1 1 1 1 0 0; ...
                         0 0 0 1 1 1 1 1 0 0 0;];
            mask = (circConv2(single(mask), conv_mask) > 0); % enlarge a bit the mask

            % Get DTM tropospheric correction
            dtm(dtm < 0) = 0; % do not consider bathimetry
            if nargout < 5 && nargin < 3
                h_correction = Core_Utils.interp1LS([h_o; 5000 * ones(100,1)], [med_tropo;  zeros(100,1)], degree, 0);
            else
                h_list = 0 : max(ceil(dtm(:)));
                h_correction = Core_Utils.interp1LS([h_o; 5000 * ones(100,1)], [med_tropo;  zeros(100,1)], degree, h_list);
            
                % Compute map of tropo corrections for height displacements
                tropo_height_correction = nan(size(dtm), 'single');
                tropo_height_correction(:) = single(h_correction(round(max(0, dtm(:)) + 1)) - h_correction(1));
            end
            
            % List of valide epochs (opening an aproximate window around the points)
            x_list = x_mg(mask);
            y_list = y_mg(mask);
            
            epoch = 1 : s_time.length();
            if numel(rate) == 1
                epoch_list = 1 : ss_rate : numel(epoch);
            else
                epoch_list = rate;
            end
            
            tropo_grid = nan(size(mask,1), size(mask,2), numel(epoch_list), 'single'); 
            
            if flag_show
                % IMAGE DEBUG: 
                fig_handle = figure; maximizeFig(fig_handle);
            end

            w_bar = Core.getWaitBar();
            w_bar.createNewBar('Generating maps of troposphere');
            w_bar.setBarLen(numel(epoch_list));
            
            tropo_clim = tropo_res + h_correction(1);
            tropo_clim = [perc(tropo_clim(:),0.0025) perc(tropo_clim(:),0.9975)];
            tropo_clim(2,:) = [perc(tropo(:),0.0025) perc(tropo(:),0.9975)];

            if flag_show
                % subplot(1,2,1);
                imh = imagesc(x_grid, y_grid, tropo_height_correction);
                if FTP_Downloader.checkNet()
                    plot_google_map('alpha', 0.65, 'MapType', 'satellite');
                end
                xlabel('Longitude [deg]');
                ylabel('Latitude [deg]');
                caxis(tropo_clim(1,:));
                %cmap = Cmap.get('c51',512);
                %colormap(flipud(cmap(2:end,:)));
                %colormap(Cmap.noaaRain);
                colormap(flipud(Cmap.get('RdBu')));
                colorbar;
                th = title(sprintf([par_name ' [cm] map @%s at sea level'], time.getEpoch(1).toString('yyyy-mm-dd HH:MM:SS')), 'FontSize', 22);
                    
                % ax2 = subplot(1,2,2);
                % imh2 = imagesc(x_grid, y_grid, tropo_height_correction);
                % if FTP_Downloader.checkNet()
                %     plot_google_map('alpha', 0.65, 'MapType', 'satellite');
                %     %plot_google_map('alpha', 0.65, 'MapType', 'roadmap');
                % end
                % xlabel('Longitude [deg]');
                % ylabel('Latitude [deg]');
                % caxis(tropo_clim(2,:));
                % cmap = Cmap.get('c51', 501);
                % colormap(flipud(cmap(2:end,:)));
                % colorbar;
                % th2 = title(ax2, 'at ground level', 'FontSize', 22);                
            end

            for i = 1 : numel(epoch_list)
                e = epoch_list(i);
                if sum(id_ok(:, epoch(e))) > 2
                    th.String = sprintf([par_name ' [cm] map %s at sea level'], time.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS'));
                    tmp = nan(size(mask));
                    if strmatch(method, 'fun')
                        tmp(mask) = funInterp2(x_list, y_list, xyu(id_ok(:, epoch(e)),1), xyu(id_ok(:, epoch(e)),2), tropo_res(epoch(e), id_ok(:, epoch(e)))', fun);
                    else
                        finterp = scatteredInterpolant(xyu(id_ok(:, epoch(e)),1),xyu(id_ok(:, epoch(e)),2), tropo_res(epoch(e), id_ok(:, epoch(e)))', method, 'none');
                        tmp(mask) = finterp(x_list, y_list);
                    end
                    tropo_grid(:,:,i) = single(tmp) + h_correction(1);
                    if flag_show                        
                        imh.CData = tropo_grid(:,:,i);
                        imh.AlphaData = ~isnan(tropo_grid(:,:,i));
                        %imh2.CData = tropo_grid(:,:,i) + tropo_height_correction;
                        %imh2.AlphaData = ~isnan(tropo_grid(:,:,i));
                        drawnow;
                    end
                else
                    if i > 1
                        tropo_grid(:,:,i) = tropo_grid(:,:,i) - 1;
                    end
                end
                w_bar.goTime(i);
            end
            w_bar.close();
        end
        
        function [tropo_out, tropo_height_correction, time] = getTropoInterp(sta_list, par_name, dlat_out, dlon_out, h_out, rate)
            % Get interpolated map of tropospheric parameter
            % Resolution is determined by the dtm in use (2 * 0.029 degrees)
            % The map is computer only on ground (> 10m) + 1 degree of margin
            %
            % INPUT
            %   par_name    type of tropospheric parameter:
            %                - ztd
            %                - zwd
            %                - pwv
            %   dlat_out    latitude of the point of interpolation  [deg]
            %   dlon_out    longitude of the point of interpolation [deg]
            %   h_out       orthometric height of the point of interpolation [m]            
            %   rate        rate in seconds, nearest to closest observation 
            %               it should be a subsample of the data rate (e.g. 300 with 30s data)
            %
            % SYNTAX
            %   [tropo_out, tropo_height_correction, time] = sta_list.getTropoInterp(par_name, dlat_out, dlon_out, h_out, rate)
            
            % Defining interpolation
            if numel(sta_list) > 2
                method = 'natural';
            else
                method = 'linear';
            end
            fun = @(dist) 0.2 * exp(-(dist)*1e1);% + 0*exp(-(dist*5e1).^2);
            %fun = @(dist) 1./(dist+1e-5);
            
            sta_list = sta_list(~sta_list.isEmptyOut_mr);
            switch lower(par_name)
                case 'ztd'
                    [tropo, s_time] = sta_list.getZtd_mr();
                case 'zwd'
                    [tropo, s_time] = sta_list.getZwd_mr();
                case 'dztd'
                    [tropo, s_time] = sta_list.getZtd_mr();
                    s_time = s_time.getEpoch(2:s_time.length);
                    s_time = s_time.addIntSeconds(s_time.getRate / 2);
                    tropo = (diff(tropo, 1)*1e2) / s_time.getRate * 3600;
                case 'dzwd'
                    [tropo, s_time] = sta_list.getZwd_mr();
                    s_time = s_time.getEpoch(2:s_time.length);
                    s_time = s_time.addIntSeconds(s_time.getRate / 2);
                    tropo = (diff(tropo, 1)*1e2) / s_time.getRate * 3600;
                case 'gn'
                    [~, s_time, ~, ~, tropo] = sta_list.getZwd_mr();
                case 'ge'
                    [~, s_time, ~, tropo, ~]  = sta_list.getZwd_mr();
                case 'dir'
                    [~, s_time, ~, ge, gn]  = sta_list.getZwd_mr();
                    tropo = atan2d(gn, ge) - 90;
                case 'pwv'
                    [tropo, s_time] = sta_list.getPwv_mr();
                case 'zhd'
                case 'nsat'
            end
            tropo = tropo * 1e2;
            
            if nargin == 6
                t_ref = round(s_time.getMatlabTime * 86400 / s_time.getRate) * s_time.getRate;
                start = t_ref(1);
                stop = t_ref(end);
                start = round(start / rate) * rate;
                stop = round(stop / rate) * rate;
                time = start : rate : stop;
                [t_ref, id_subset] = intersect(t_ref, time);
                time = GPS_Time(t_ref / 86400, [], true);
            else
                time = s_time;
                id_subset =  1 : time.length();
            end
            med_tropo = mean(tropo, 1, 'omitnan')';

            degree = 4;
            coo = Coordinates.fromXYZ(sta_list.getMedianPosXYZ);
            [lat, lon, up, h_o] = coo.getGeodetic;
            
            xyu = [dlon_out, dlat_out, h_out];
            
            % Remove missing stations
            id_ok = ~isnan(tropo'); % logical matrix of valid tropo
                                                
            % Correct data for height
            
            h_correction_ztd = sta_list.getZtdReduction(2, [], false) * 1e2;
            tropo_res = bsxfun(@minus, tropo', h_correction_ztd)';            
            
            tropo_height_correction = sta_list.getZtdReduction(2, [(dlon_out * (pi/180)) .* cos((dlat_out * (pi/180))), (dlat_out * (pi/180)) h_out], false) * 1e2;
                        
            % List of valide epochs (opening an aproximate window around the points)
            x_list = lon ./ pi * 180;
            y_list = lat ./ pi * 180;
            
            epoch = 1 : time.length();
            epoch_list = 1 : numel(epoch);
            
            tropo_out = nan(size(dlat_out,1), numel(epoch_list), 'single'); 
                     
            w_bar = Core.getWaitBar();
            w_bar.createNewBar('Generating maps of troposphere');
            w_bar.setBarLen(numel(epoch_list));
            
            for i = 1 : numel(epoch_list)
                e = epoch_list(i);
                if sum(id_ok(:, epoch(e))) > 1
                    if strcmp(method, 'fun')
                        tmp = funInterp2(dlon_out, dlat_out, x_list(id_ok(:, id_subset(epoch(e)))), y_list(id_ok(:, id_subset(epoch(e)))), tropo_res(id_subset(epoch(e)), id_ok(:, id_subset(epoch(e))))', fun);
                    else
                        warning off
                        finterp = scatteredInterpolant(x_list(id_ok(:, id_subset(epoch(e)))), y_list(id_ok(:, id_subset(epoch(e)))), tropo_res(id_subset(epoch(e)), id_ok(:, id_subset(epoch(e))))', method, 'linear');
                        warning on
                        tmp = finterp(dlon_out, dlat_out);
                    end
                    try
                        tropo_out(:,i) = single(tmp);
                    catch % faster then if checking if tmp is empty
                        tmp = funInterp2(dlon_out, dlat_out, x_list(id_ok(:, id_subset(epoch(e)))), y_list(id_ok(:, id_subset(epoch(e)))), tropo_res(id_subset(epoch(e)), id_ok(:, id_subset(epoch(e))))', fun);
                        tropo_out(:,i) = tmp;
                    end
                else
                    if sum(id_ok(:, id_subset(epoch(e)))) == 1
                        tropo_out(:,i) = single(tropo_res(id_subset(epoch(e)), id_ok(:, id_subset(epoch(e)))));
                    else
                        tropo_out(:,i) = single(nan);
                        % If I have no station return nan
                    end
                end
                w_bar.goTime(i);
            end
            w_bar.close();
            tropo_out = tropo_out';
        end

        function [tropo_out, tropo_height_correction, time] = getTropoInterpOld(sta_list, par_name, dlat_out, dlon_out, h_out)
            % Get interpolated map of tropospheric parameter
            % Resolution is determined by the dtm in use (2 * 0.029 degrees)
            % The map is computer only on ground (> 10m) + 1 degree of margin
            %
            % INPUT
            %   par_name    type of tropospheric parameter:
            %                - ztd
            %                - zwd
            %                - pwv
            %   rate        rate in seconds, nearest to closest observation 
            %               it should be a subsample of the data rate (e.g. 300 with 30s data)
            %   flag_show   if true show debug images
            %
            % SYNTAX
            %   [tropo_grid, x_grid, y_grid, time, tropo_height_correction] = sta_list.getTropoMap(par_name, rate)
            
            % Defining interpolation
            method = 'natural';
            fun = @(dist) 0.2 * exp(-(dist)*1e1) + 0*exp(-(dist*5e1).^2);
            %fun = @(dist) 1./(dist+1e-5);
            
            sta_list = sta_list(~sta_list.isEmptyOut_mr);
            switch lower(par_name)
                case 'ztd'
                    [tropo, s_time] = sta_list.getZtd_mr();
                case 'zwd'
                    [tropo, s_time] = sta_list.getZwd_mr();
                case 'dztd'
                    [tropo, s_time] = sta_list.getZtd_mr();
                    s_time = s_time.getEpoch(2:s_time.length);
                    s_time = s_time.addIntSeconds(s_time.getRate / 2);
                    tropo = (diff(tropo, 1)*1e2) / s_time.getRate * 3600;
                case 'dzwd'
                    [tropo, s_time] = sta_list.getZwd_mr();
                    s_time = s_time.getEpoch(2:s_time.length);
                    s_time = s_time.addIntSeconds(s_time.getRate / 2);
                    tropo = (diff(tropo, 1)*1e2) / s_time.getRate * 3600;
                case 'gn'
                    [~, s_time, ~, ~, tropo] = sta_list.getZwd_mr();
                case 'ge'
                    [~, s_time, ~, tropo, ~]  = sta_list.getZwd_mr();
                case 'dir'
                    [~, s_time, ~, ge, gn]  = sta_list.getZwd_mr();
                    tropo = atan2d(gn, ge) - 90;
                case 'pwv'
                    [tropo, s_time] = sta_list.getPwv_mr();
                case 'zhd'
                case 'nsat'
            end
            tropo = tropo * 1e2;
            
            time = s_time; % s_time.getEpoch(1 : s_time.length());
            med_tropo = mean(tropo, 1, 'omitnan')';

            degree = 4;
            coo = Coordinates.fromXYZ(sta_list.getMedianPosXYZ);
            [lat, lon, up, h_o] = coo.getGeodetic;
            
            xyu = [dlon_out, dlat_out, h_out];
            
            % Remove missing stations
            id_ok = ~isnan(tropo'); % logical matrix of valid tropo
                                                
            % Correct data for height
            tropo_height = Core_Utils.interp1LS(h_o, med_tropo, degree);
            tropo_res = bsxfun(@minus, tropo', tropo_height)';            
            
            h_correction = Core_Utils.interp1LS([h_o; 5000 * ones(100,1)], [med_tropo;  zeros(100,1)], degree, [0; h_out(:)]);
            
            % Compute map of tropo corrections for height displacements
            tropo_height_correction = h_correction(2:end) - h_correction(1);
                        
            % List of valide epochs (opening an aproximate window around the points)
            x_list = lon ./ pi * 180;
            y_list = lat ./ pi * 180;
            
            epoch = 1 : s_time.length();
            epoch_list = 1 : numel(epoch);
            
            tropo_out = nan(size(dlat_out,1), numel(epoch_list), 'single'); 
                     
            w_bar = Core.getWaitBar();
            w_bar.createNewBar('Generating maps of troposphere');
            w_bar.setBarLen(numel(epoch_list));
            
            for i = 1 : numel(epoch_list)
                e = epoch_list(i);
                if sum(id_ok(:, epoch(e))) > 2
                    if strmatch(method, 'fun')
                        tmp = funInterp2(dlon_out, dlat_out, x_list(id_ok(:, epoch(e))), y_list(id_ok(:, epoch(e))), tropo_res(epoch(e), id_ok(:, epoch(e)))', fun);
                    else
                        finterp = scatteredInterpolant(x_list(id_ok(:, epoch(e))), y_list(id_ok(:, epoch(e))), tropo_res(epoch(e), id_ok(:, epoch(e)))', method, 'linear');
                        tmp = finterp(dlon_out, dlat_out);
                    end
                    tropo_out(:,i) = single(tmp) + h_correction(1);
                else
                    if i > 1
                        tropo_out(:,i) = tropo_out(:,i) - 1;
                    end
                end
                w_bar.goTime(i);
            end
            w_bar.close();
            tropo_out = tropo_out';
        end

        function [id_rec, dist_3d, dist_up] = getCloserRec(sta_list, lat, lon, h_o)
            % Get the id and 3D distance of the closest station w.r.t. a given point
            %
            % INPUT 
            %   lat, lon    [degree]
            %   h_o         [m] orthometric height
            %
            % OUTPUT
            %   id_rec      id of the closest GNSS station
            %   dist_3d     minimum distance 3D [m] w.r.t. the requested coordinates
            %   dist_up     minimum distance up [m] w.r.t. the requested coordinates
            %
            % SYNTAX
            %   [id_rec, dist_3d, dist_up] = sta_list.getCloserRec(lat, lon, h_o)
            
            if nargin == 3
                h_o = 0;
            end
            
            sta_xyz = sta_list.getMedianPosXYZ;
            out_coo = Coordinates.fromGeodetic(lat / 180 * pi, lon / 180 * pi, [], h_o);
            out_xyz = out_coo.getXYZ;
            
            n_out = size(out_xyz, 1);
            n_sta = size(sta_xyz, 1);
            
            % check all the distances
            d2 = (repmat(sta_xyz(:,1), 1, n_out) - repmat(out_xyz(:,1)', n_sta, 1)) .^2 + ...
                (repmat(sta_xyz(:,2), 1, n_out) - repmat(out_xyz(:,2)', n_sta, 1)) .^2 + ...
                (repmat(sta_xyz(:,3), 1, n_out) - repmat(out_xyz(:,3)', n_sta, 1)) .^2;
            
            % Keep the closest station
            [d2_min, id_rec] = min(d2);
            dist_3d = sqrt(d2_min);
            [~, ~, ~, h_ortho] = sta_list(id_rec).getMedianPosGeodetic();
            dist_up = h_ortho - h_o;
        end
        
        function rec_works = getWork(sta_list, id)
            % Return the working receiver for a GNSS station array
            %
            % SYNTAX
            %  rec_works = sta_list.getWork(<id>)
            if nargin < 2
                rec_works = [sta_list.work];
            else
                id(id > numel(sta_list)) = [];
                rec_works = [sta_list(id).work];
            end
        end
        
        function [id_rds, lat, lon, appr_dist] = getCloseRaobIdList(sta_list, max_dist)
            % Get all the radiosondes ids and their positions
            % that are located close to the GNSS stations
            %
            % INPUT
            %   max_dist      [ km ]
            %
            % SYNTAX
            %   [id_rds, lat_lon_rds] = sta_list.getCloseRadiosondeId(max_dist)
            id_rds = cell(numel(sta_list), 1);
            [lat, lon, appr_dist] = deal(nan(numel(sta_list), 1));
            for s = 1 : numel(sta_list)
                [lat_tmp, lon_tmp] = sta_list(s).getMedianPosGeodetic();
                [id_rds(s), lat(s), lon(s)] = Radiosonde.getCloseStations(lat_tmp, lon_tmp, 1);
                coo_gnss = Coordinates.fromGeodetic(lat_tmp./180*pi, lon_tmp./180*pi, 0);
                coo_raob =  Coordinates.fromGeodetic(lat(s)./180*pi, lon(s)./180*pi, 0);
                raob_dist = coo_gnss - coo_raob;
                appr_dist(s) = sqrt(sum(raob_dist.xyz_diff.^2, 2)) .* 1e-3; % distance in km
            end
            % Sort them by distance, it's a trick to keep only
            % the closer distances
            [appr_dist, id_sort] = sort(appr_dist);
            id_rds = id_rds(id_sort);            
            [id_rds, id_ok] = unique(id_rds);
            appr_dist = appr_dist(id_ok);
            lat = lat(id_sort(id_ok));
            lon = lon(id_sort(id_ok));  
            if nargin == 2
                id_ok = appr_dist <= max_dist;
                appr_dist = appr_dist(id_ok);
                id_rds = id_rds(id_ok);
                lat = lat (id_ok);
                lon = lon(id_ok);
            end
        end
    end

    % ==================================================================================================================================================
    %% SETTER
    % ==================================================================================================================================================
    methods (Access = public)
        function marker_name = setMarkerName(this, marker_name)
            % Set the Marker name of the Receiver
            %
            % SYNTAX:
            %   marker_name = setMarkerName(this)
            this.marker_name = marker_name;
            if ~isempty(this.work)
                this.work.rinex_file_name = marker_name(1 : max(4, numel(marker_name))); % Force 4Ch name to match this
            end                
        end
        
        function setSlantFilterSize(this, win_size)
            % Setter multi_receiver to change filtering windows size for slant export
            %
            % SYNTAX
            %   this.setSlantFilterSize(win_size)
            for r = 1 : numel(this)
                this(r).slant_filter_win = win_size;
            end
        end
        
        function fixPos(sta_list, mode)
            % fix the position of the receiver into the reference frame
            % object
            %
            % SYNTAX:
            %  sta_list.fixPos(mode)
            for s = 1 : length(sta_list)
                if strcmpi(mode,'work') % get from work
                    xyz = sta_list(s).work.getMedianPosXYZ();
                    Core.getReferenceFrame.setCoo(sta_list(s).getMarkerName4Ch, xyz, 2, [0 0 0], [0 0], GPS_Time([1970 1 1 0 0 0]), GPS_Time([2099 1 1 0 0 0]));
                elseif strcmpi(mode,'out')
                    xyz = sta_list(s).out.getPosXYZ();
                    xyz = xyz(end,:);
                    Core.getReferenceFrame.setCoo(sta_list(s).getMarkerName4Ch, xyz, 2, [0 0 0], [0 0], GPS_Time([1970 1 1 0 0 0]), GPS_Time([2099 1 1 0 0 0]));
                end
            end
        end
        
        function unFixPos(sta_list)
            % unfix the position of the receiver into the reference frame
            % object
            %
            % SYNTAX:
            %  sta_list.unFixPos()
            for s = 1 : length(sta_list)
                Core.getReferenceFrame.setFlag(sta_list(s).getMarkerName4Ch, 1)
            end
        end
    end

    % ==================================================================================================================================================
    %% TESTER
    % ==================================================================================================================================================
    methods (Access = public)
        function checkTropo_mr(sta_list, thr_lev)
            % Check ZWD for possible outliers (works on the mr ZWD getter)
            % And remove all the outliers from the receiver tropo parameters 
            %
            % INPUT
            %   thr_lev       level of threshold (usually 5.5 - conservative)
            %
            % SYNTAX
            %   id_ko = sta_list.checkTropo_mr(thr_lev, flag_reduce)
            
            % state = Core.getCurrentSettings; 
            if nargin < 2 || isempty(thr_lev)
                thr_lev = [];
            end
            
            log = Core.getLogger();
            if (isempty(sta_list))
                log.addWarning('No valid receiver found');
            else
                log.addMarkedMessage('Cleaning tropospheric parameters');
                out_list = [sta_list.out];
                out_list = out_list(~sta_list.isEmptyOut_mr);
                if isempty(out_list)
                    work_list = [sta_list.work];
                    work_list = work_list(~sta_list.isEmptyWork_mr);
                    work_list.cleanTropo(thr_lev)
                else
                    out_list.cleanTropo(thr_lev);
                end
                log.addStatusOk('Cleaning done');
            end
        end
    end
    
    % ==================================================================================================================================================
    %% STATIC FUNCTIONS used as utilities
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function [p_time, id_sync] = getSyncTimeExpanded(rec, p_rate, use_pos_time)
            % Get the common time among all the receivers
            %
            % SYNTAX
            %   [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(rec, <p_rate>, <use_pos_time>);
            %
            % EXAMPLE:
            %   [p_time, id_sync] = GNSS_Station.getSyncTimeExpanded(rec, 30);

            if nargin < 3 || isempty(use_pos_time)
                use_pos_time = false;
            end

            if sum(~rec.isEmpty_mr) == 0
                % no valid receiver
                p_time = GPS_Time;
                id_sync = [];
            else
                if nargin < 2 || isempty(p_rate)
                    p_rate = 1e-6;

                    for r = 1 : numel(rec)
                        if (rec(r).out.time.length) > 2
                            if use_pos_time
                                p_rate = lcm(round(p_rate * 1e6), round(rec(r).out.time_pos.getRate * 1e6)) * 1e-6; % enable this line to sync rates
                            else
                                p_rate = lcm(round(p_rate * 1e6), round(rec(r).out.time.getRate * 1e6)) * 1e-6; % enable this line to sync rates
                            end
                        end
                    end
                end

                % prepare reference time
                % processing time will start with the receiver with the last first epoch
                %          and it will stop  with the receiver with the first last epoch

                out = [rec.out];
                first_id_ok = find(~out.isEmpty_mr, 1, 'first');
                if ~isempty(first_id_ok)
                    if use_pos_time
                        p_time_zero = round(rec(first_id_ok).out.time_pos.first.getMatlabTime() * 24)/24; % get the reference time for positions
                    else
                        p_time_zero = round(rec(first_id_ok).out.time.first.getMatlabTime() * 24)/24; % get the reference time
                    end
                end

                % Get all the common epochs
                t = [];
                for r = 1 : numel(rec)
                    if use_pos_time
                        rec_rate = min(86400, iif(rec(r).out.time_pos.length == 1, 86400, rec(r).out.time_pos.getRate));
                        t = [t; round(rec(r).out.time_pos.getRefTime(p_time_zero) / rec_rate) * rec_rate];
                    else
                        rec_rate = min(5, rec(r).out.time.getRate);
                        t = [t; round(rec(r).out.time.getRefTime(p_time_zero) / rec_rate) * rec_rate];
                    end
                    % p_rate = lcm(round(p_rate * 1e6), round(rec(r).out.time.getRate * 1e6)) * 1e-6; % enable this line to sync rates
                end
                t = unique(t);

                % If p_rate is specified use it
                if nargin > 1
                    t = intersect(t, (t(1) : p_rate : t(end) + p_rate)');
                end

                % Create reference time
                p_time = GPS_Time(p_time_zero, t);
                id_sync = nan(p_time.length(), numel(rec));

                % Get intersected times
                for r = 1 : numel(rec)
                    if use_pos_time
                        rec_rate = iif(rec(r).out.time_pos.length == 1, 86400, rec(r).out.time_pos.getRate);
                        [~, id1, id2] = intersect(t, round(rec(r).out.time_pos.getRefTime(p_time_zero) / rec_rate) * rec_rate);
                    else
                        rec_rate = min(5, rec(r).out.time.getRate);
                        [~, id1, id2] = intersect(t, round(rec(r).out.time.getRefTime(p_time_zero) / rec_rate) * rec_rate);
                    end

                    id_sync(id1, r) = id2;
                end
            end
        end

        function [p_time, id_sync] = getSyncTimeTR(sta_list, obs_type, p_rate)
            % Get the common (shortest) time among all the used receivers and the target(s)
            % For each target (obs_type == 0) produce a different cell array with the sync of the other receiver
            % e.g.  Reference receivers @ 1Hz, trg1 @1s trg2 @30s
            %       OUTPUT 1 sync @1Hz + 1 sync@30s
            %
            % SYNTAX
            %   [p_time, id_sync] = Receiver.getSyncTimeTR(rec, obs_type, <p_rate>);
            %
            % SEE ALSO:
            %   this.getSyncTimeExpanded
            %
            if nargin < 3
                p_rate = 1e-6;
            end
            if nargin < 2
                % choose the longest as reference
                len = zeros(1, numel(sta_list));
                for r = 1 : numel(sta_list)
                    len(r) = sta_list(r).out.length;
                end
                obs_type = ones(1, numel(sta_list));
                obs_type(find(len == max(len), 1, 'first')) = 0;
            end

            % Do the target(s) as last
            [~, id] = sort(obs_type, 'descend');

            % prepare reference time
            % processing time will start with the receiver with the last first epoch
            %          and it will stop  with the receiver with the first last epoch

            first_id_ok = find(~sta_list.isEmptyOut_mr, 1, 'first');
            p_time_zero = round(sta_list(first_id_ok).out.time.first.getMatlabTime() * 24)/24; % get the reference time
            p_time_start = sta_list(first_id_ok).out.time.first.getRefTime(p_time_zero);
            p_time_stop = sta_list(first_id_ok).out.time.last.getRefTime(p_time_zero);
            p_rate = lcm(round(p_rate * 1e6), round(sta_list(first_id_ok).out.time.getRate * 1e6)) * 1e-6;

            p_time = GPS_Time(); % empty initialization

            i = 0;
            for r = id
                ref_t{r} = sta_list(r).out.time.getRefTime(p_time_zero);
                if obs_type(r) > 0 % if it's not a target
                    if ~sta_list(r).out.isEmpty
                        p_time_start = max(p_time_start,  round(sta_list(r).out.time.first.getRefTime(p_time_zero) * sta_list(r).out.time.getRate) / sta_list(r).out.time.getRate);
                        p_time_stop = min(p_time_stop,  round(sta_list(r).out.time.last.getRefTime(p_time_zero) * sta_list(r).out.time.getRate) / sta_list(r).out.time.getRate);
                        p_rate = lcm(round(p_rate * 1e6), round(sta_list(r).out.time.getRate * 1e6)) * 1e-6;
                    end
                else
                    % It's a target

                    % recompute the parameters for the ref_time estimation
                    % not that in principle I can have up to num_trg_rec ref_time
                    % in case of multiple targets the reference times should be independent
                    % so here I keep the temporary rt0 rt1 r_rate var
                    % instead of ref_time_start, ref_time_stop, ref_rate
                    pt0 = max(p_time_start, round(sta_list(r).out.time.first.getRefTime(p_time_zero) * sta_list(r).out.time.getRate) / sta_list(r).out.time.getRate);
                    pt1 = min(p_time_stop, round(sta_list(r).out.time.last.getRefTime(p_time_zero) * sta_list(r).out.time.getRate) / sta_list(r).out.time.getRate);
                    pr = lcm(round(p_rate * 1e6), round(sta_list(r).out.time.getRate * 1e6)) * 1e-6;
                    pt0 = ceil(pt0 / pr) * pr;
                    pt1 = floor(pt1 / pr) * pr;

                    % return one p_time for each target
                    i = i + 1;
                    p_time(i) = GPS_Time(p_time_zero, (pt0 : pr : pt1));
                    p_time(i).toUnixTime();

                    id_sync{i} = nan(p_time(i).length, numel(id));
                    for rs = id % for each rec to sync
                        if ~sta_list(rs).out.isEmpty && ~(obs_type(rs) == 0 && (rs ~= r)) % if it's not another different target
                            [~, id_ref, id_rec] = intersect(round(sta_list(rs).out.time.getRefTime(p_time_zero) * 1e5)/1e5, (pt0 : pr : pt1));
                            id_sync{i}(id_rec, rs) = id_ref;
                        end
                    end
                end
            end
        end

        function bsl_ids = getBaselineId(n_rec)
            % Get id of all the combinations of the stations
            %
            % SYNTAX:
            %   bsl_id = GNSS_Station.getBaselineId(n_rec);
            [r1, r2] = meshgrid(1 : n_rec, 1 : n_rec);
            bsl_ids = [serialize(tril(r1, -1)) serialize(tril(r2, -1))];
            bsl_ids = bsl_ids(bsl_ids(:, 1) > 0 & bsl_ids(:, 2) > 0, :);
        end
        
        function ant_mp = getCurrentMPM(applied_ant, type)
            % Get the currently applied multipath map 
            %
            % OUTPUT
            %   mp.type              type of applied multipath (if 0 no mp)
            %   mp.time_lim          reference time of the applied map
            %   mp.(sys_c).(trk).map current map according to Core.getState.flag_rec_mp
            %
            % SYNTAX 
            %   [applied_map, type] = getAppliedMP(this, ant_model, <type>)
            
            if nargin < 2 || isempty(type)
                ant_mp.type = Core.getState.flag_rec_mp; % get the kind of applied map
            else
                ant_mp.type = type;
            end
            
            try
                ant_mp.time_lim = applied_ant.time_lim;
                % Get the satellite systems available in the zerniche multipath struct
                sys_c_list = intersect(cell2mat(fields(applied_ant)'), 'GRECJI');
                for sys_c = sys_c_list
                    if isfield(applied_ant, sys_c)
                        % This constellation is already present into the applied Zernike MultiPath
                        trk_list = fields(applied_ant.(sys_c))';
                        for trk = trk_list
                            trk = trk{1};
                            if ~isfield(ant_mp, sys_c)
                                ant_mp.(sys_c) = struct;
                            end
                            if ~isfield(ant_mp.(sys_c), trk)
                                ant_mp.(sys_c).(trk) = struct;
                            end
                            switch ant_mp.type
                                case 1, ant_mp.(sys_c).(trk).map = applied_ant.(sys_c).(trk).z_map;
                                case 2, ant_mp.(sys_c).(trk).map = applied_ant.(sys_c).(trk).r_map;
                                case 3, ant_mp.(sys_c).(trk).map = applied_ant.(sys_c).(trk).g_map;
                                case 4, ant_mp.(sys_c).(trk).map = applied_ant.(sys_c).(trk).c_map;
                                case 5, ant_mp.(sys_c).(trk).map = applied_ant.(sys_c).(trk).g1_map;
                                case 6, ant_mp.(sys_c).(trk).map = applied_ant.(sys_c).(trk).c1_map;
                            end
                        end
                    end
                end
            catch ex
                if ant_mp.type == 0
                   % this could be normal 
                else
                    ant_mp.type = 0;
                end
            end
        end
    end
    %% STATIC FUNCTIONS used as tools
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function showGlobalMap(new_fig, flag_labels, flag_polar)
            % Show Map of the stations
            % downloading the DTM and showing it
            %
            % CITATION:
            %   Pawlowicz, R., 2019. "M_Map: A mapping package for MATLAB", version 1.4k, [Computer software],
            %   available online at www.eoas.ubc.ca/~rich/map.html.
            %
            % INPUT
            %   new_fig     open a new figure
            %   flag_labels show/not show labels
            %   flag_polar  could be 'N' / 'S' / false
            %
            % SYNTAX
            %   GNSS_Station.showGlobalMap(new_fig, flag_labels, flag_polar);
            %
            % EXAMPLE
            %   GNSS_Station.showGlobalMap()
            
            if nargin < 2 || isempty(flag_labels)
                flag_labels = false;
            end
            if nargin < 3 || isempty(flag_polar)
                flag_polar = false;
            end
            flag_large_points = true;
            point_size = 20;
            large_point_size = iif(flag_polar, 50, 25);
            point_color = [14, 25, 41]/256;
            resolution = 'high';
            
            if nargin < 2 || isempty(new_fig)
                new_fig = true;
            end
            if new_fig
                f = figure('Visible', 'off');
            else
                f = gcf;
                hold on;
            end
            
            fh_list = f;
            fig_name = sprintf('WorldRecMapDtm');
            f.UserData = struct('fig_name', fig_name);
                                    
            Core.getLogger.addMarkedMessage('Preparing map, please wait...');
            
            maximizeFig(f);
            f.Color = [1 1 1];
            
            
            [coo, name] = GNSS_Station.getStationFromWorldArchive();
            [lat, lon, h] =  coo.getGeodetic;
            lat = lat / pi * 180;
            lon = lon / pi * 180;
                        
            % set map limits
            if numel(name) == 1
                lon_lim = minMax(lon) + [-0.05 0.05];
                lat_lim = minMax(lat) + [-0.05 0.05];
            else
                lon_lim = minMax(lon); lon_lim = lon_lim + [-1 1] * diff(lon_lim)/15;
                lat_lim = minMax(lat); lat_lim = lat_lim + [-1 1] * diff(lat_lim)/15;
            end
            
            lon_lim = max(-179.999, min(179.999, lon_lim));
            lat_lim = max(-89.999, min(89.999, lat_lim));
            
            lon_lim = lon_lim + [-1 1] * max(0, (0.5 - diff(minMax(lon_lim))) / 2);
            lat_lim = lat_lim + [-1 1] * max(0, (0.5 - diff(minMax(lat_lim))) / 2);
            
            nwse = [lat_lim(2), lon_lim(1), lat_lim(1), lon_lim(2)];
            clon = nwse([2 4]) + [-0.001 0.001];
            clat = max(-90, min(90, nwse([3 1]) + [-0.001 0.001]));
            clon = max(-180, min(180, clon));
            
            if (flag_polar)
                if flag_polar == 'S'                    
                    m_proj('stereographic','lat',-90,'long',0,'radius', 25);
                    id_ko = lat > -65;
                else
                    m_proj('stereographic','lat',90,'long',0,'radius', 25);
                    id_ko = lat < 65;
                end
                lat(id_ko) = [];
                lon(id_ko) = [];
                name = name(~id_ko);
            else
                %m_proj('equidistant','lon',clon,'lat',clat);   % Projection
                m_proj('miller', 'lon', clon, 'lat', clat);   % Projection
            end
            
            axes
            cmap = flipud(gray(1000)); colormap(cmap(150: end, :));
            
            % retrieve external DTM
            if flag_polar
                % use ETOPO instead
                colormap(gca, cmap(100 : end - 100, :))
                m_etopo2('shadedrelief', 'gradient', 3);
            else
                try
                    cmap = flipud(gray(1000)); colormap(gca, cmap(150: end, :));
                    [dtm, lat_dtm, lon_dtm] = Core.getRefDTM(nwse, 'ortho', resolution);

                    sensor_pole = sum(~isnan(dtm) / size(dtm, 2),2);
                    id_pole = (find(sensor_pole > 0.5, 1, 'first') - 1);
                    if any(id_pole)
                        dtm(1 : id_pole, :) = repmat(dtm(id_pole + 1, : ), id_pole, 1);
                    end
                    dtm = flipud(dtm);
                    id_pole = (find(flipud(sensor_pole) > 0.5, 1, 'first') - 1);
                    if any(id_pole)
                        dtm(1 : id_pole, :) = repmat(dtm(id_pole + 1, : ), id_pole, 1);
                    end
                    
                    % comment the following line to have bathimetry
                    dtm(dtm < 0) = nan; % - 1/3 * max(dtm(:));
                    
                    % uncomment the following line to have colors
                    %colormap(gca, Cmap.adaptiveTerrain(minMax(dtm(:))));
                    drawnow;
                    
                    warning off;
                    [shaded_dtm, x, y] = m_shadedrelief(lon_dtm, lat_dtm, dtm, 'nan', [0.98, 0.98, 1], 'gra', 30);
                    warning on;
                    %h_dtm = m_pcolor(lon_dtm, lat_dtm, dtm);
                    %h_dtm.CData = shaded_dtm;
                    
                    m_image(lon_dtm, lat_dtm, shaded_dtm);
                catch ex
                    % use ETOPO1 instead
                    colormap(gca, cmap(100 : end, :))
                    m_etopo2('shadedrelief','gradient', 3);
                end
            end
            
            % read shapefile
            shape = 'none';
            if (~strcmp(shape,'none'))
                if (~strcmp(shape,'coast')) && (~strcmp(shape,'fill'))
                    if (strcmp(shape,'10m'))
                        M = m_shaperead('countries_10m');
                    elseif (strcmp(shape,'30m'))
                        M = m_shaperead('countries_30m');
                    else
                        M = m_shaperead('countries_50m');
                    end
                    [x_min, y_min] = m_ll2xy(min(lon_lim), min(lat_lim));
                    [x_max, y_max] = m_ll2xy(max(lon_lim), max(lat_lim));
                    for k = 1 : length(M.ncst)
                        lam_c = M.ncst{k}(:,1);
                        ids = lam_c <  min(lon);
                        lam_c(ids) = lam_c(ids) + 360;
                        phi_c = M.ncst{k}(:,2);
                        [x, y] = m_ll2xy(lam_c, phi_c);
                        if sum(~isnan(x))>1
                            x(find(abs(diff(x)) >= abs(x_max - x_min) * 0.90) + 1) = nan; % Remove lines that occupy more than th 90% of the plot
                            line(x,y,'color', [0.3 0.3 0.3]);
                        end
                    end
                else
                    if (strcmp(shape,'coast'))
                        m_coast('line','color', lineCol);
                    else
                        m_coast('patch',lineCol);
                    end
                end
            end
            
            hold on;
            
            if (flag_polar)
                if flag_polar == 'S'
                    m_grid('XaxisLocation', 'top', 'tickdir','in', 'fontsize', 16);
                else
                    m_grid('tickdir','in', 'fontsize', 16);
                end
                % m_ruler(1.1, [.05 .40], 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
                ylim([-0.47 0.47]);
                drawnow
            else
                m_grid('box','fancy','tickdir','in', 'fontsize', 16);
                % m_ruler(1.1, [.05 .40], 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
                drawnow
                % m_ruler([.7 1], -0.05, 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
            end
            
            [x, y] = m_ll2xy(lon, lat);
            %point_color = Cmap.get('viridis', numel(x));
            %point_size = 25;
            if size(point_color, 1) > 1
                scatter(x(:), y(:), point_size, 1:numel(x), 'filled'); hold on;
                colormap(point_color);
            else
                plot(x(:), y(:),'.', 'MarkerSize', point_size, 'Color', point_color); hold on;
            end
            if flag_labels
                % Label BG (in background w.r.t. the point)
                for r = 1 : numel(x)
                    tmp = name{r};
                    txt = text(x(r), y(r), ['     ' name ' '], ...
                        'FontWeight', 'bold', 'FontSize', 12, 'Color', [1 1 1], ...
                        'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                        'Margin', 2, 'LineWidth', 2, ...
                        'HorizontalAlignment','left');
                end
            end
            
            if flag_large_points
                plot(x(:), y(:), 'ko', 'MarkerSize', large_point_size/3, 'LineWidth', 1 + 1 * logical(flag_polar));
                clim = minMax(h(:));
                n_col = ceil(diff(clim)/10);
                color = round((h(:) - clim(1)) / 10) + 1;
                colormap(Cmap.linspaced);
                cb = colorbar; title(cb, '[m]');
                caxis(clim);                
                for r = 1 : numel(x)                    
                    plot(x(r), y(r), '.', 'MarkerSize', large_point_size, 'Color', Core_UI.getColor(color(r), n_col));
                end
                plot(x(:), y(:), '.k', 'MarkerSize', large_point_size/7);
            end
            
            if flag_labels
                for r = 1 : numel(x)
                    tmp = name{r};
                    t = text(x(r), y(r), ['   ' tmp], ...
                        'FontWeight', 'bold', 'FontSize', 12, 'Color', [0 0 0], ...
                        ...%'FontWeight', 'bold', 'FontSize', 10, 'Color', [0 0 0], ...
                        ...%'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                        'Margin', 2, 'LineWidth', 2, ...
                        'HorizontalAlignment','left');
                end
            end
            Core_UI.addExportMenu(f); Core_UI.addBeautifyMenu(f); Core_UI.beautifyFig(f);
            f.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
            title(sprintf('Stations position\\fontsize{5} \n'), 'FontSize', 16);
            %xlabel('Longitude [deg]');
            %ylabel('Latitude [deg]');
            ax = gca; ax.FontSize = 16;
            Core.getLogger.addStatusOk('The map is ready ^_^');
        end
        
        function [coo, marker_name, flag] = getStationFromWorldArchive(marker_filter)
            % Read the coordinates present in:
            %   INSTALL_DIR/../data/station/CRD/gred_world_archive.xyz
            %
            % INPUT 
            %   marker_filter      when it is specified a filter is applied
            %
            % SYNTAX
            %   [coo, marker_name, flag] = GNSS_Station.getStationFromWorldArchive(<marker_name>);
            
            crd_dir = File_Name_Processor.getFullDirPath(fullfile(Core.getInstallDir(), '../data/station/CRD/'));
            file_name = fullfile(crd_dir, 'gred_world_archive.xyz');
            if ~exist(file_name, 'file')
                Core.getLogger.addError(sprintf('"%s" not found', file_name));
                crd_list = [];
            else
                [txt, lim] = Core_Utils.readTextFile(file_name, 3);
                
                if isempty(txt)
                    crd_list = [];
                else
                    name_list = txt(repmat(lim(3:end,1),1,19) + repmat((0:18),size(lim,1)-2,1));
                    flag = txt(lim(3:end,1) + 61);
                    
                    % strip unique coordinates from marker names
                    id_sta = flag == '1';
                    name_list(id_sta, :) = reshape(regexprep(serialize(name_list(id_sta, :)')', '\_\d{4}\_\d{4}', '          '), 19, sum(id_sta))';
                    
                    % coordinates
                    id_tmp = repmat(lim(3 : end,1),1,42) + repmat((19:60),size(lim,1)-2,1);
                    xyz = reshape(sscanf(serialize(txt(id_tmp)')', '%f'), 3, size(lim, 1) -2)';
                    
                    % remove recievers with 0 coordinate
                    id_sta = flag == '0';
                    name_list(id_sta, :) = [];
                    xyz(id_sta, :) = [];
                    flag(id_sta) = [];
                                        
                    marker_name = {};
                    id_ok = true(size(name_list, 1), 1);
                    for s = 1 : size(name_list, 1)
                        marker_name{s} = strtrim(name_list(s, :));
                        if nargin == 1 && ~isempty(marker_filter)
                            id_ok(s) = contains(marker_name{s}, marker_filter);
                        end
                    end
                    
                    marker_name = marker_name(id_ok);
                    xyz = xyz(id_ok, :);
                    flag = flag(id_ok);
                    coo = Coordinates.fromXYZ(xyz);
                    
                    [~, ~, h] =  coo.getGeodetic;
                    id_ok = (h > -200); % Coordinates with elevation < 200m are definitely not valid
                    
                    if any(~id_ok)
                        marker_name = marker_name(id_ok);
                        xyz = xyz(id_ok, :);
                        flag = flag(id_ok);
                        coo = Coordinates.getElement(id_ok);
                    end
                end
            end
        end
        
        function [coo, marker_list, flag, d] = getCloseStations(lat, lon, n_stations)
            % Get the marker name of the stations in a certain area,
            % or close to a point
            %
            % SYNTAX
            %   [coo, marker_list, flag, d] = GNSS_Station.getCloseStations(lat_lim, lon_lim);
            %   [coo, marker_list, flag, d] = GNSS_Station.getCloseStations(lat, lon, n_stations)
            %   [coo, marker_list, flag, d] = GNSS_Station.getCloseStations(world_marker, [], n_stations)
            %
            % EXAMPLE
            %   [coo, marker_list, flag, d] = GNSS_Station.getCloseStations([45.25 45.75], [9 9.5]);
            %   [coo, marker_list, flag, d] = GNSS_Station.getCloseStations(45.43, 9.28, 3)
            %   [coo, marker_list, flag, d] = GNSS_Station.getCloseStations('MIL2', [], 3)
            
            if nargin >= 1 && ischar(lat)
                % in this case I search for coordinates in the World map
                marker_filter = lat;
                [coo, marker_name] = GNSS_Station.getStationFromWorldArchive(marker_filter);
                
                if isempty(marker_name)
                    Core.getLogger.addWarning(sprintf('No station containing "%s" found', marker_filter))
                else
                    if numel(marker_name) > 1
                        Core.getLogger.addWarning(sprintf('Multiple stations containing "%s" found\nChosing the first\n Check GNSS_Station.getCloseStations()', marker_filter));
                    end
                    [lat, lon] =  coo.getElement(1).getGeodetic();
                    lat = lat / pi * 180;
                    lon = lon / pi * 180;
                end
            end
            
            [coo, marker_list, flag] = GNSS_Station.getStationFromWorldArchive();
            [lat_sta, lon_sta, h] =  coo.getGeodetic;
            lat_sta = lat_sta / pi * 180;
            lon_sta = lon_sta / pi * 180;
            
            if numel(lat) > 1
                % lat is a limit
                lat = sort(lat);
                lon = sort(lon);

                id_ok = (lat_sta >= lat(1)) & (lat_sta <= lat(2)) & ...
                    (lon_sta >= lon(1)) & (lon_sta <= lon(2));
                if (nargout == 4)
                    d = sphericalDistance(lat_sta, lon_sta, lat, lon);
                    d = d(id_ok);
                end
                marker_list = marker_list(id_ok);
                coo = coo.getElement(id_ok);
            else
                % distance from a coordinate
                d = sphericalDistance(lat_sta, lon_sta, lat, lon);
                [d, id_ok] = sort(d);
                marker_list = marker_list(id_ok);
                coo = coo.getElement(id_ok);
                flag = flag(id_ok);
            end
            
            if nargin == 3 && ~isempty(n_stations)
                id_ok = 1 : min(n_stations, length(marker_list));
                d = d(id_ok);
                marker_list = marker_list(id_ok);
                coo = coo.getElement(id_ok);
                flag = flag(id_ok);                
            end            
        end
    end
    %% METHODS PLOTTING FUNCTIONS
    % ==================================================================================================================================================

    % Various debug images
    % name variant:
    %   c cartesian
    %   s scatter
    %   p polar
    %   m mixed
    methods (Access = public)
        function showAll(sta_list)
            % Try to show all the possible plots            
            for i = 1:numel(sta_list)
                sta_list(i).out.showAll;
            end
        end

        function fh_list = showObsStats(sta_list)
            % Show statistics about the observations stored in the object
            %
            % SYNTAX
            %   this.showObsStats()

            fh_list = [];
            for s = 1 : numel(sta_list)
                fh_list = [fh_list; sta_list(s).work.showObsStats()]; %#ok<AGROW>
            end
        end

        function fh_list = showDataAvailability(sta_list, sys_list)
            % Show data availability for each receiver workspace
            %
            % SYNTAX
            %   this.showDataAvailability(sys_list)
            
            sta_list = sta_list(not(sta_list.isEmpty_mr));
            fh_list = [];
            if nargin < 2
                sys_list = Core.getConstellationCollector.getAvailableSys();
            end
            for s = 1 : numel(sta_list)
                fh_list = [fh_list; sta_list(s).work.showDataAvailability(sys_list)];
            end
        end
        
        function fh_list = showProcessingQualityInfo(sta_list)
            % Show quality info about the receiver processing
            % SYNTAX this.showProcessingQualityInfo();

            fh_list = [];
            for r = 1 : length(sta_list)
                if ~sta_list(r).isEmpty() && ~sta_list(r).out.isEmpty()
                    rec_out = sta_list(r).out;
                    fh_list = [fh_list; rec_out.showProcessingQualityInfo()]; %#ok<AGROW>
                end
            end
        end

        function fh_list = showSyncronizedPos(sta_list)
            % Show syncronized 
            rec = sta_list(not(sta_list.isEmpty_mr));
            rate = rec(1).getPos.time.getRate;
            t0 = inf;
            tend = 0;
            for r = 1 : numel(rec)
                coo = rec(r).getPos;
                t = floor(coo.time.getMatlabTime()*(86400/rate));
                tc{r} = t;
                
                t0 = min(t0, tc{r}(1));
                tend = max(tend, tc{r}(1));
            end
            t_sync = nan(tend - t0 + 1, numel(rec));
            xyz_sync = nan(tend - t0 + 1, 3, numel(rec));
            for r = 1 : numel(rec)
                t_sync(tc{r} - t0 + 1, r) = r + (mod(r, 2)*2)-1;
                xyz_sync(tc{r} - t0 + 1, :, r) = rec(r).getPos.getXYZ;
            end
            
            % realign XYZ coordinates
            for c = 1 : 3
                tmp = squeeze(zero2nan(xyz_sync(:,c,:)));
                [~, id_max] = max(sum(not(isnan(tmp))));
                med_bsl = median(zero2nan(tmp) - tmp(:, id_max), 'omitnan');
                tmp = squeeze(tmp - med_bsl);
                tmp = bsxfun(@minus, tmp, median(tmp, 2, 'omitnan'));
                tmp = tmp + med_bsl;
                xyz_sync(:,c,:) = reshape(tmp, size(tmp,1), 1, size(tmp,2));
            end
            
            % XYZ to ENU
            enu_sync = [];
            for r = 1: numel(rec)
                coo = Coordinates.fromXYZ(xyz_sync(:,:,r));
                enu_sync(:,:,r) = coo.getLocal(coo.getMedianPos);
            end
            
            fh = figure;
            %imagesc(((1:size(t_sync,1)) + t0) /(86400/rate) , 1:size(t_sync,2), t_sync'); colormap(Cmap.get('viridis',20));
            tmp_e = squeeze(zero2nan(enu_sync(:,1,:)));
            tmp_n = squeeze(zero2nan(enu_sync(:,2,:)));
            tmp_u = squeeze(zero2nan(enu_sync(:,3,:)));
            ax_list = [];
            
            subplot(3,1,1); imagesc(1:size(t_sync,1), 1:size(t_sync,2), sqrt(tmp_e.^2)');
            title('Sessions with a position');
            ax = gca;
            ax.YAxis.TickValues = 1 : numel(rec);
            for r = 1 : numel(rec)
                ax.YAxis.TickLabels{r} = rec(r).getMarkerName4Ch;
            end
            ax_list(1) = ax;
            caxis([-0.03 0.05])
            ylabel('East')
            
            subplot(3,1,2); imagesc(1:size(t_sync,1), 1:size(t_sync,2), sqrt(tmp_n.^2)');
            ax = gca;
            ax.YAxis.TickValues = 1 : numel(rec);
            for r = 1 : numel(rec)
                ax.YAxis.TickLabels{r} = rec(r).getMarkerName4Ch;
            end
            ax_list(2) = ax;
            caxis([-0.03 0.05])
            ylabel('North')
            
            subplot(3,1,3); imagesc(1:size(t_sync,1), 1:size(t_sync,2), sqrt(tmp_u.^2)');
            ax = gca;
            ax.YAxis.TickValues = 1 : numel(rec);
            for r = 1 : numel(rec)
                ax.YAxis.TickLabels{r} = rec(r).getMarkerName4Ch;
            end
            ax_list(3) = ax;
            caxis([-0.03 0.05])
            
            linkaxes(ax_list);
            %imagesc(1:size(t_sync,1), 1:size(t_sync,2), sqrt(tmp_x.^2 + tmp_y.^2 + tmp_z.^ 2)');
            colormap(Cmap.get('viridis',20));
            xlabel('Session')
            ylabel('Up')
            
            fig_name = sprintf('Coo_report');
            fh.UserData = struct('fig_name', fig_name);
            Core_UI.beautifyFig(fh);
            Core_UI.addExportMenu(fh);
            Core_UI.addBeautifyMenu(fh);
            fh_list = fh;
        end
        
        function fh_list = showPositionENU(sta_list, flag_add_coo, n_obs)
            % Plot East North Up coordinates of the receiver
            %
            % INPUT:
            %   sta_list                 list of GNSS_Station objects
            %   flag_add_coo             use external / internal / high rate cordinates
            %   n_obs                    use only the last n_obs
            %
            % SYNTAX 
            %   this.plotPositionENU(<flag_add_coo>, <n_obs>);

            if nargin < 3 || isempty(n_obs) || isnan(n_obs)
                n_obs = 0;
            end
            if ~(nargin >= 2 && ~isempty(flag_add_coo) && flag_add_coo ~= 0)
                flag_add_coo = 0;
            end
            
            if flag_add_coo >= 0
                sta_list = sta_list(~sta_list.isEmptyOut_mr);
            end
            out_list = [sta_list.out];
            if numel(out_list) > 0
                fh_list = out_list.showPositionENU(flag_add_coo, n_obs);
            end
        end

        function fh_list = showPositionPlanarUp(sta_list, flag_add_coo, n_obs)
            % Plot Planar and Up coordinates of the receiver
            %
            % INPUT:
            %   sta_list                 list of GNSS_Station objects
            %   flag_add_coo             use external / internal / high rate cordinates
            %   n_obs                    use only the last n_obs
            %
            % SYNTAX 
            %   this.showPositionPlanarUp(<flag_add_coo>, <n_obs>);
            
            if nargin < 3 || isempty(n_obs) || isnan(n_obs)
                n_obs = 0;
            end
            if ~(nargin >= 2 && ~isempty(flag_add_coo) && flag_add_coo ~= 0)
                flag_add_coo = 0;
            end
            
            if flag_add_coo >= 0
                sta_list = sta_list(~sta_list.isEmptyOut_mr);
            end
            out_list = [sta_list.out];
            fh_list = out_list.showPositionPlanarUp(flag_add_coo, n_obs);
        end
        
        function fh_list = showPositionXYZ(sta_list, flag_add_coo, n_obs)
            % Plot X Y Z coordinates of the receiver
            %
            % INPUT:
            %   sta_list                 list of GNSS_Station objects
            %   flag_add_coo             use external / internal / high rate cordinates
            %   n_obs                    use only the last n_obs
            %
            % SYNTAX this.plotPositionXYZ(<flag_add_coo>, <n_obs>);
            if nargin < 3 || isempty(n_obs) || isnan(n_obs)
                n_obs = 0;
            end
            if ~(nargin >= 2 && ~isempty(flag_add_coo) && flag_add_coo ~= 0)
                flag_add_coo = 0;
            end
            
            if flag_add_coo >= 0
                sta_list = sta_list(~sta_list.isEmptyOut_mr);
            end
            out_list = [sta_list.out];
            fh_list = out_list.showPositionXYZ(flag_add_coo, n_obs);
        end

        function fh_list = showPositionSigmas(sta_list)
            % Show Sigmas of the solutions
            %
            % SYNTAX
            %   this.showPositionSigmas();

            fh_list = [];
            for r = 1 : length(sta_list)
                rec = sta_list(r).out;
                if ~isempty(rec)
                    fh_list = [fh_list; rec.showPositionSigmas()]; %#ok<AGROW>
                end
            end
        end

        function fh_list = showMap(sta_list, new_fig)
            % Show Google Map of the stations
            %
            % CITATION:
            %   Pawlowicz, R., 2019. "M_Map: A mapping package for MATLAB", version 1.4k, [Computer software],
            %   available online at www.eoas.ubc.ca/~rich/map.html.
            %
            % INPUT
            %   new_fig     open a new figure
            %
            % SYNTAX
            %   sta_list.showMap(new_fig);
            if nargin < 2
                new_fig = true;
            end
            fh_list = sta_list.showMapGoogle(new_fig);
        end

        function fh_list = showMapDtm(sta_list, new_fig, resolution, add_lat, add_lon)
            % Show Map of the stations
            % downloading the DTM and showing it
            %
            % CITATION:
            %   Pawlowicz, R., 2019. "M_Map: A mapping package for MATLAB", version 1.4k, [Computer software],
            %   available online at www.eoas.ubc.ca/~rich/map.html.
            %
            % INPUT
            %   new_fig     open a new figure
            %   resolution  'high' / 'low'
            %
            % SYNTAX
            %   sta_list.showMapDtm(new_fig, resolution);
            
            N_MAX_STA = 50;
            sta_list = sta_list(~sta_list.isEmpty_mr);
            
            flag_labels = true;
            flag_large_points = true;
            point_size = 15;
            %point_color = [14, 25, 41]/256;
            point_color = Core_UI.ORANGE;
            
            if nargin < 2 || isempty(new_fig)
                new_fig = true;
            end
            if new_fig
                f = figure('Visible', 'on');
            else
                f = gcf;
                hold on;
            end
            
            fh_list = f;
            fig_name = sprintf('RecMapDtm');
            f.UserData = struct('fig_name', fig_name);
            
            if (nargin < 3) || isempty(resolution)
                resolution = 'low';
            end
            % check accepted values (low / high)
            switch resolution
                case 'high'
                otherwise
                    resolution = 'low';
            end
            Core.getLogger.addMarkedMessage('Preparing map, please wait...');
            
            f.Color = [1 1 1];
            [lat, lon] = sta_list.getMedianPosGeodetic();

            if nargin == 5
                lat_tmp = [lat(:); add_lat(:)];
                lon_tmp = [lon(:); add_lon(:)];
            else
                lat_tmp = lat;
                lon_tmp = lon;
            end
            % set map limits
            if numel(lon_tmp) == 1
                lon_lim = minMax(lon_tmp) + [-0.05 0.05];
                lat_lim = minMax(lat_tmp) + [-0.05 0.05];
            else
                lon_lim = minMax(lon_tmp); lon_lim = min(179.999, max(-179.999, lon_lim + [-1 1] * diff(lon_lim) / 6));
                lat_lim = minMax(lat_tmp); lat_lim = min(89.999, max(-89.999,lat_lim + [-1 1] * diff(lat_lim) / 6));
            end
            lon_lim = lon_lim + [-1 1] * max(0, (0.5 - diff(minMax(lon_lim))) / 2);
            lat_lim = lat_lim + [-1 1] * max(0, (0.5 - diff(minMax(lat_lim))) / 2);
            nwse = [lat_lim(2), lon_lim(1), lat_lim(1), lon_lim(2)];
            clon =  min(180, max(-180, nwse([2 4]) + [-0.001 0.001]));
            clat = max(-90, min(90, nwse([3 1]) + [-0.001 0.001]));

            m_proj('equidistant','lon',clon,'lat',clat);   % Projection
            %m_proj('utm', 'lon',lon_lim,'lat',lat_lim);   % Projection
            drawnow; axes
            cmap = flipud(gray(1000)); colormap(cmap(150: end, :));
            
            % retrieve external DTM
            try
                [dtm, lat_dtm, lon_dtm] = Core.getRefDTM(nwse, 'ortho', resolution);
                dtm = flipud(dtm);
                % comment the following line to have bathimetry
                dtm(dtm < 0) = nan; % - 1/3 * max(dtm(:));
            
                % uncomment the following line to have colors
                % colormap(Cmap.adaptiveTerrain(minMax(dtm(:))));
                drawnow;
                
                [shaded_dtm, x, y] = m_shadedrelief(lon_dtm, lat_dtm, dtm, 'nan', [0.98, 0.98, 1], 'gra', 30);
                %h_dtm = m_pcolor(lon_dtm, lat_dtm, dtm);
                %h_dtm.CData = shaded_dtm;
                m_image(lon_dtm, lat_dtm, shaded_dtm); 
            catch
                % use ETOPO1 instead
                m_etopo2('shadedrelief','gradient', 3);
            end

            % read shapefile
            shape = 'none';
            if (~strcmp(shape,'none'))
                if (~strcmp(shape,'coast')) && (~strcmp(shape,'fill'))
                    if (strcmp(shape,'10m'))
                        M = m_shaperead('countries_10m');
                    elseif (strcmp(shape,'30m'))
                        M = m_shaperead('countries_30m');
                    else
                        M = m_shaperead('countries_50m');
                    end
                    [x_min, y_min] = m_ll2xy(min(lon_lim), min(lat_lim));
                    [x_max, y_max] = m_ll2xy(max(lon_lim), max(lat_lim));
                    for k = 1 : length(M.ncst)
                        lam_c = M.ncst{k}(:,1);
                        ids = lam_c <  min(lon);
                        lam_c(ids) = lam_c(ids) + 360;
                        phi_c = M.ncst{k}(:,2);
                        [x, y] = m_ll2xy(lam_c, phi_c);
                        if sum(~isnan(x))>1
                            x(find(abs(diff(x)) >= abs(x_max - x_min) * 0.90) + 1) = nan; % Remove lines that occupy more than th 90% of the plot
                            line(x,y,'color', [0.3 0.3 0.3]);
                        end
                    end
                else
                    if (strcmp(shape,'coast'))
                        m_coast('line','color', lineCol);
                    else
                        m_coast('patch',lineCol);
                    end
                end
            end
            
            hold on;
            
            m_grid('box','fancy','tickdir','in', 'fontsize', Core_UI.getFontSize(11));
            % m_ruler(1.1, [.05 .40], 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
            drawnow
            m_ruler([.7 1], -0.075, 'tickdir','out','ticklen',[.007 .007], 'fontsize',Core_UI.getFontSize(9));
            [x, y] = m_ll2xy(lon, lat);
            
            %point_color = Cmap.get('viridis', numel(x));
            %point_size = 25;
            if size(point_color, 1) > 1
                scatter(x(:), y(:), point_size, 1:numel(x), 'filled'); hold on;
                colormap(point_color);
            else
                if numel(sta_list) > N_MAX_STA
                    plot(x(:), y(:),'.', 'MarkerSize', point_size, 'Color', Core_UI.BLACK); hold on;
                    plot(x(:), y(:),'.', 'MarkerSize', point_size-3, 'Color', point_color); hold on;
                else
                    plot(x(:), y(:),'.', 'MarkerSize', 8, 'Color', Core_UI.BLACK); hold on;
                end
            end
            if flag_labels
                % Label BG (in background w.r.t. the point)
                if numel(sta_list) < N_MAX_STA
                    for r = 1 : numel(sta_list)
                        name = upper(sta_list(r).getMarkerName4Ch());
                        text(x(r), y(r), ['     ' name ' '], ...
                            'FontWeight', 'bold', 'FontSize', 12, 'Color', [1 1 1], ...
                            'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                            'Margin', 2, 'LineWidth', 2, ...
                            'HorizontalAlignment','left');
                    end
                end
            end
            
            if flag_large_points && numel(sta_list) < N_MAX_STA
                for r = 1 : numel(sta_list)
                    plot(x(r), y(r), '.', 'MarkerSize', 45, 'Color', Core_UI.getColor(r, numel(sta_list)), 'UserData', 'GNSS_point');
                end
                plot(x(:), y(:), '.k', 'MarkerSize', 5);
                plot(x(:), y(:), 'ko', 'MarkerSize', 15, 'LineWidth', 2);
            end
            
            if flag_labels
                if numel(sta_list) < N_MAX_STA
                    for r = 1 : numel(sta_list)
                        name = upper(sta_list(r).getMarkerName4Ch());
                        text(x(r), y(r), ['   ' name], ...
                            'FontWeight', 'bold', 'FontSize', 12, 'Color', [0 0 0], ...
                            'Margin', 2, 'LineWidth', 2, ...
                            'HorizontalAlignment','left');
                    end
                end
            end
            Core_UI.addExportMenu(f); Core_UI.addBeautifyMenu(f); Core_UI.beautifyFig(f, 'light');
            f.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
            title(sprintf('Map of GNSS stations\\fontsize{5} \n'), 'FontSize', 16);
            %xlabel('Longitude [deg]');
            %ylabel('Latitude [deg]');
            ax = gca; ax.FontSize = 16;
            Core.getLogger.addStatusOk('The map is ready ^_^');
        end

        function fh_list = showMapGoogle(sta_list, new_fig, add_lat, add_lon)
            % Show Google Map of the stations
            %
            % CITATION:
            %   Pawlowicz, R., 2019. "M_Map: A mapping package for MATLAB", version 1.4k, [Computer software],
            %   available online at www.eoas.ubc.ca/~rich/map.html.
            %
            % INPUT
            %   new_fig     open a new figure
            %
            % SYNTAX
            %   sta_list.showMapGoogle(new_fig);
           
            N_MAX_STA = 50;
            sta_list = sta_list(~sta_list.isEmpty_mr);
            
            flag_labels = true;
            flag_large_points = true;
            point_size = 15;
            point_color = [0, 255, 10]/256;
            
            if nargin < 2 || isempty(new_fig)
                new_fig = true;
            end
            if new_fig
                f = figure('Visible', 'off');
            else
                f = gcf;
                hold on;
            end
            
            fh_list = f;
            fig_name = sprintf('RecMapGoogle');
            f.UserData = struct('fig_name', fig_name);

            Core.getLogger.addMarkedMessage('Preparing map, please wait...');
            
            f.Color = [1 1 1];
            [lat, lon] = sta_list.getMedianPosGeodetic();
            
            if nargin == 4
                lat_tmp = [lat(:); add_lat(:)];
                lon_tmp = [lon(:); add_lon(:)];
            else
                lat_tmp = lat;
                lon_tmp = lon;
            end
            % set map limits
            if numel(lon_tmp) == 1
                lon_lim = minMax(lon_tmp) + [-0.05 0.05];
                lat_lim = minMax(lat_tmp) + [-0.05 0.05];
            else
                lon_lim = minMax(lon_tmp); lon_lim = lon_lim + [-1 1] * diff(lon_lim) / 6;
                lat_lim = minMax(lat_tmp); lat_lim = lat_lim + [-1 1] * diff(lat_lim) / 6;
            end
            lon_lim = lon_lim + [-1 1] * max(0, (0.5 - diff(minMax(lon_lim))) / 2);
            lat_lim = lat_lim + [-1 1] * max(0, (0.5 - diff(minMax(lat_lim))) / 2);
            nwse = [lat_lim(2), lon_lim(1), lat_lim(1), lon_lim(2)];
            clon = nwse([2 4]) + [-0.001 0.001];
            clat = max(-90, min(90, nwse([3 1]) + [-0.001 0.001]));

            axes
            xlim(clon);
            ylim(clat);
            [lon_ggl,lat_ggl, img_ggl] = plot_google_map('alpha', 0.95, 'maptype','satellite','refresh',0,'autoaxis',0);
            xlim(lon_lim);
            ylim(lat_lim);
            
            m_proj('equidistant','lon',clon,'lat',clat);   % Projection
            %m_proj('utm', 'lon',lon_lim,'lat',lat_lim);   % Projection
            drawnow
            m_image(lon_ggl, lat_ggl, img_ggl);
            
            % read shapefile
            shape = 'none';
            if (~strcmp(shape,'none'))
                if (~strcmp(shape,'coast')) && (~strcmp(shape,'fill'))
                    if (strcmp(shape,'10m'))
                        M = m_shaperead('countries_10m');
                    elseif (strcmp(shape,'30m'))
                        M = m_shaperead('countries_30m');
                    else
                        M = m_shaperead('countries_50m');
                    end
                    [x_min, y_min] = m_ll2xy(min(lon_lim), min(lat_lim));
                    [x_max, y_max] = m_ll2xy(max(lon_lim), max(lat_lim));
                    for k = 1 : length(M.ncst)
                        lam_c = M.ncst{k}(:,1);
                        ids = lam_c <  min(lon);
                        lam_c(ids) = lam_c(ids) + 360;
                        phi_c = M.ncst{k}(:,2);
                        [x, y] = m_ll2xy(lam_c, phi_c);
                        if sum(~isnan(x))>1
                            x(find(abs(diff(x)) >= abs(x_max - x_min) * 0.90) + 1) = nan; % Remove lines that occupy more than th 90% of the plot
                            line(x,y,'color', [0.3 0.3 0.3]);
                        end
                    end
                else
                    if (strcmp(shape,'coast'))
                        m_coast('line','color', lineCol);
                    else
                        m_coast('patch',lineCol);
                    end
                end
            end
            hold on;
            
            m_grid('box','fancy','tickdir','in', 'fontsize', 16);
            % m_ruler(1.1, [.05 .40], 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
            drawnow
            %m_ruler([.7 1], -0.05, 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
            [x, y] = m_ll2xy(lon, lat);
            
            
            %point_color = Cmap.get('viridis', numel(x));
            %point_size = 25;
            if size(point_color, 1) > 1
                scatter(x(:), y(:), point_size, 1:numel(x), 'filled'); hold on;
                colormap(point_color);
            else
                if numel(sta_list) > N_MAX_STA
                    plot(x(:), y(:),'.', 'MarkerSize', point_size, 'Color', Core_UI.BLACK); hold on;
                    plot(x(:), y(:),'.', 'MarkerSize', point_size-3, 'Color', point_color); hold on;
                else
                    plot(x(:), y(:),'.', 'MarkerSize', 8, 'Color', Core_UI.BLACK); hold on;
                end
            end
            if flag_labels
                % Label BG (in background w.r.t. the point)
                if numel(sta_list) < N_MAX_STA
                    for r = 1 : numel(sta_list)
                        name = upper(sta_list(r).getMarkerName4Ch());
                        txt = text(x(r), y(r), ['   ' name], ...
                            'FontWeight', 'bold', 'FontSize', 12, 'Color', [1 1 1], ...
                            'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                            'Margin', 2, 'LineWidth', 2, ...
                            'HorizontalAlignment','left');
                    end
                end
            end
            
            if flag_large_points && numel(sta_list) < N_MAX_STA
                for r = 1 : numel(sta_list)
                    plot(x(r), y(r), '.', 'MarkerSize', 45, 'Color', Core_UI.getColor(r, numel(sta_list)), 'UserData', 'GNSS_point');
                end
                plot(x(:), y(:), '.k', 'MarkerSize', 5);
                plot(x(:), y(:), 'ko', 'MarkerSize', 15, 'LineWidth', 2);
            end
            if flag_labels
                if numel(sta_list) < N_MAX_STA
                    for r = 1 : numel(sta_list)
                        name = upper(sta_list(r).getMarkerName4Ch());
                        t = text(x(r), y(r), ['   ' name], ...
                            'FontWeight', 'bold', 'FontSize', 12, 'Color', [0 0 0], ...
                            ...%'FontWeight', 'bold', 'FontSize', 10, 'Color', [0 0 0], ...
                            ...%'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                            'Margin', 2, 'LineWidth', 2, ...
                            'HorizontalAlignment','left');
                        %t.Units = 'pixels';
                        %t.Position(1) = t.Position(1) + 20 + 10 * double(numel(sta_list) == 1);
                        %t.Units = 'data';
                    end
                end
            end
            
            Core_UI.addExportMenu(f); Core_UI.addBeautifyMenu(f); Core_UI.beautifyFig(f);
            f.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
            title(sprintf('Map of GNSS stations\\fontsize{5} \n'), 'FontSize', 16);
            %xlabel('Longitude [deg]');
            %ylabel('Latitude [deg]');
            ax = gca; ax.FontSize = 16;
            Core.getLogger.addStatusOk('The map is ready ^_^');
        end
        
        function fh_list = showMapGoogleLegacy(sta_list, new_fig)
            % Show Google Map of the stations
            % Old version without m_map
            %
            % SYNTAX
            %   sta_list.showMapGoogle(new_fig);
            if nargin < 2
                new_fig = true;
            end
            if new_fig
                fh = figure('Visible', 'off');
            else
                fh = gcf;
                hold on;
            end
            
            fh_list = fh;
            fig_name = sprintf('RecMapLgc');
            fh.UserData = struct('fig_name', fig_name);

            [lat, lon] = sta_list.getMedianPosGeodetic();

            plot(lon(:), lat(:),'.k', 'MarkerSize', 5); hold on;            
            % Label BG (in background w.r.t. the point)
            for r = 1 : numel(sta_list)
                name = upper(sta_list(r).getMarkerName4Ch());
                text(lon(r), lat(r), ['   ' name], ...
                    'FontWeight', 'bold', 'FontSize', 12, 'Color', [1 1 1], ...
                    'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                    'Margin', 2, 'LineWidth', 2, ...
                    'HorizontalAlignment','left');
            end
            
            for r = 1 : numel(sta_list)
                plot(lon(r), lat(r), '.', 'MarkerSize', 45, 'Color', Core_UI.getColor(r, numel(sta_list)), 'UserData', 'GNSS_point');
            end
            plot(lon(:), lat(:), '.k', 'MarkerSize', 5);
            plot(lon(:), lat(:), 'ko', 'MarkerSize', 15, 'LineWidth', 2);

            if numel(sta_list) == 1
                lon_lim = minMax(lon);
                lat_lim = minMax(lat);
                lon_lim(1) = lon_lim(1) - 0.05;
                lon_lim(2) = lon_lim(2) + 0.05;
                lat_lim(1) = lat_lim(1) - 0.05;
                lat_lim(2) = lat_lim(2) + 0.05;
            else
                lon_lim = xlim();
                lon_lim(1) = lon_lim(1) - diff(lon_lim)/3;
                lon_lim(2) = lon_lim(2) + diff(lon_lim)/3;
                lat_lim = ylim();
                lat_lim(1) = lat_lim(1) - diff(lat_lim)/3;
                lat_lim(2) = lat_lim(2) + diff(lat_lim)/3;
            end
            lon_lim = lon_lim + [-1 1] * max(0, (0.5 - diff(minMax(lon_lim))) / 2);
            lat_lim = lat_lim + [-1 1] * max(0, (0.5 - diff(minMax(lat_lim))) / 2);
            
            xlim(lon_lim);
            ylim(lat_lim);

            for r = 1 : numel(sta_list)
                name = upper(sta_list(r).getMarkerName4Ch());
                t = text(lon(r), lat(r), ['   ' name], ...
                    'FontWeight', 'bold', 'FontSize', 12, 'Color', [0.1 0.1 0.1], ...
                    ...%'FontWeight', 'bold', 'FontSize', 10, 'Color', [0 0 0], ...
                    ...%'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                    'Margin', 2, 'LineWidth', 2, ...
                    'HorizontalAlignment','left');
                %t.Units = 'pixels';
                %t.Position(1) = t.Position(1) + 20 + 10 * double(numel(sta_list) == 1);
                %t.Units = 'data';
            end

            plot_google_map('alpha', 0.95, 'MapType', 'satellite');
            title(sprintf('Map of GNSS stations\\fontsize{5} \n'), 'FontSize', 16);
            xlabel('Longitude [deg]');
            ylabel('Latitude [deg]');
            Core_UI.addExportMenu(fh); Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh);
            fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
            Core.getLogger.addStatusOk('The map is ready ^_^');
        end

        function fh_list = showCMLRadarMapAniRec(sta_list)
            fh_list = [];
            try
                cml = CML();
                fh_list = cml.showRadarMapAniRec(sta_list);
            catch
                sta_list(1).log.addError('You need GReD utilities to have this features');
            end
        end
        
        function fh_list = showAniMapTropoScatter(sta_list, par_name, new_fig, epoch, flag_export)
            % Show a tropo map with all the station in sta_list
            %
            % INPUT
            %   tropo_par   accepted tropo parameter:
            %               - 'zwd'
            %               - 'ztd'
            %   new_fig     if true or missing open a new figure
            %   epoch       list of epoch to display (pay attention that there is a subsampling rate in the function 1:10:end)
            %   flag_export if true try to export a video (the frames are not going to be seen)
            %               the video file will be saved in the out folder specified in the project
            %
            % SYNTAX
            %   sta_list.showMapTropo(<new_fig>, <epoch>, <flag_export>);
            
            ss_rate = 10; % subsample rate for show;

            if nargin < 3 || isempty(new_fig)
                new_fig = true;
            end
            if new_fig
                f_handle = figure;
            else
                f_handle = gcf;
                hold on;
            end
            fh_list = f_handle;
            fig_name = sprintf('RecAniMap_%s', upper(par_name));
            f_handle.UserData = struct('fig_name', fig_name);

            if nargin < 5
                flag_export = false;
            end

            maximizeFig(f_handle);

            switch lower(par_name)
                case 'ztd'
                    [res_tropo, s_time] = sta_list.getZtd_mr();
                    par_str = 'ZTD';
                    par_str_short = 'ZTD';
                case 'ztd_red'
                    [res_tropo, s_time] = sta_list.getReducedZtd_mr();
                    par_str = 'reduced ZTD';
                    par_str_short = 'RedZTD';
                case 'zwd'
                    [res_tropo, s_time] = sta_list.getZwd_mr();
                    par_str = 'ZWD';
                    par_str_short = 'ZWD';
                case 'gn'
                case 'ge'
                case 'pwv'
                case 'zhd'
                case 'nsat'
            end
            res_tropo = res_tropo * 1e2;
                
            if nargin < 3 || isempty(epoch)
                epoch = 1 : s_time.length();
            end

            coo = Coordinates.fromXYZ(sta_list.getMedianPosXYZ);
            [lat, lon] = coo.getGeodetic;
            clear coo

            sh = scatter(lon(:)./pi*180, lat(:)./pi*180, 300, res_tropo(epoch(1),:)', 'filled');
            hold on;
            % plot(lon(:)./pi*180, lat(:)./pi*180,'ko','MarkerSize', 15, 'LineWidth', 2);
            caxis([min(res_tropo(:)) max(res_tropo(:))]);
            colormap(gat);
            colorbar;

            lon_lim = xlim();
            lon_lim(1) = lon_lim(1) - 0.1;
            lon_lim(2) = lon_lim(2) + 0.1;
            lat_lim = ylim();
            lat_lim(1) = lat_lim(1) - 0.1;
            lat_lim(2) = lat_lim(2) + 0.1;

            xlim(lon_lim);
            ylim(lat_lim);

            ax = f_handle.Children(end);
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            if new_fig
                if FTP_Downloader.checkNet()
                    plot_google_map('alpha', 0.95, 'MapType', 'satellite');
                end
                xlabel('Longitude [deg]');
                ylabel('Latitude [deg]');
            end
            th = title(sprintf([par_str ' variations [cm] map @%s'], s_time.getEpoch(1).toString('yyyy-mm-dd HH:MM:SS')), 'FontSize', 30);
            Core_UI.insertLogo(f_handle, 'SouthEast');

            if flag_export
                frameset = {};
                frame = getframe(f_handle);
                frameset{1} = frame(1:2:end,1:2:end,:); % subsample (1:2)
                f_handle.Visible = 'off';
                Core.getLogger.addMarkedMessage('Exporting video');
                fprintf('%5d/%5d', 0, 99999);
            end
            drawnow
            if numel(epoch) > 1
                epoch_list = (ss_rate + 1) : ss_rate : numel(epoch);
                for i = 1 : numel(epoch_list)
                    e = epoch_list(i);
                    if any(res_tropo(epoch(e),:))
                        th.String = sprintf([par_str ' variations [cm] map %s'], s_time.getEpoch(epoch(e)).toString('yyyy-mm-dd HH:MM:SS'));
                        sh.CData = res_tropo(epoch(e),:)';
                        if not(flag_export)
                            drawnow
                        end
                    end
                    if flag_export
                        fprintf('%s%5d/%5d',char(8 * ones(1,11)), i,numel(epoch_list));
                        frame = getframe(f_handle);
                        frameset{i + 1} = frame(1:2:end,1:2:end,:); % subsample (1:2)
                    end
                end
            end

            if flag_export
                fprintf('%s',char(8 * ones(1,11)));
                if ismac() || ispc()
                    % Better compression on Mac > 10.7 and Win > 7
                    video_out = VideoWriter(fullfile(Core.getState.getOutDir, ['AniMap' par_str_short '.mp4']), 'MPEG-4');
                else
                    % Linux doesn't have mp4 compression avaiable
                    video_out = VideoWriter(fullfile(Core.getState.getOutDir, ['AniMap' par_str_short '.avi']));
                end
                video_out.FrameRate = 30;
                video_out.Quality = 91;
                open(video_out);
                for i = 1 : numel(frameset)
                    writeVideo(video_out, frameset{i});
                end
                close(video_out);
                Core.getLogger.addStatusOk(sprintf('"%s" done ^_^', fullfile(Core.getState.getOutDir, video_out.Filename)));
                close(f_handle)
            end
        end
       
        function fh_list = showAniMapTropoInterp(sta_list_full, par_name, nwse, rate, flag_dtm, flag_export)
            % Show a tropo map with all the station in sta_list
            %
            % INPUT
            %   tropo_par   accepted tropo parameter:
            %               - 'zwd'
            %               - 'ztd'
            %   nwse        Nort West South East coordinates
            %   rate        rate in seconds, nearest to closest observation 
            %               it should be a subsample of the data rate (e.g. 300 with 30s data)
            %   flag_dtm    flag to add height_correction (default == false)
            %   flag_export if true try to export a video (the frames are not going to be seen)
            %               the video file will be saved in the out folder specified in the project
            %
            % SYNTAX
            %   sta_list.showAniMapTropoInterp(par_name, <new_fig>, <rate>, <flag_dtm>, <flag_export>);
            %
            % EXAMPLE
            %   % over Japan
            %   sta_list.showAniMapTropoInterp('ZWD', [45.8, 123.5, 23, 146.5], 200, 2, false);
            
            sta_list_full = sta_list_full(~sta_list_full.isEmptyOut_mr);
            
            switch lower(par_name)
                case 'ztd'
                    par_str = 'ZTD';
                    par_str_short = 'ZTD';
                case 'zwd'
                    par_str = 'ZWD';
                    par_str_short = 'ZWD';
                case 'gn'
                case 'ge'
                case 'pwv'
                    par_str = 'PWV';
                    par_str_short = 'PWV';
                case 'zhd'
                case 'nsat'
            end                        

            f = figure;
            fh_list = f;
            fig_name = sprintf('RecAniMapInterp_%s', upper(par_name));
            f.UserData = struct('fig_name', fig_name);
            
            if nargin < 4 || isempty(rate)
                rate = 300;
            end
            if nargin < 5
                flag_dtm = false;
            end
            if nargin < 6
                flag_export = false;
            end

            maximizeFig(f);
            
            f.Visible = 'off';
            f.Color = [1 1 1];

            % Set map projection / limits
            
            [lat, lon] = sta_list_full.getMedianPosGeodetic();
            if nargin >= 3 && ~isempty(nwse)
                margin = 0.5;
                id_ok = (lat >= (nwse(3) - margin)) & (lat <= (nwse(1) + margin)) & ...
                    (lon >= (nwse(2) - margin)) & (lon <= (nwse(4) + margin));
                sta_list = sta_list_full(id_ok);
            end
            
            if nargin < 3 || isempty(nwse)
                sta_list = sta_list_full;
                % set map limits
                if numel(sta_list) == 1
                    lon_lim = minMax(lon) + [-0.05 0.05];
                    lat_lim = minMax(lat) + [-0.05 0.05];
                else
                    lon_lim = minMax(lon); lon_lim = lon_lim + [-1 1] * diff(lon_lim)/15;
                    lat_lim = minMax(lat); lat_lim = lat_lim + [-1 1] * diff(lat_lim)/15;
                end
                lon_lim = lon_lim + [-1 1] * max(0, (0.5 - diff(minMax(lon_lim))) / 2);
                lat_lim = lat_lim + [-1 1] * max(0, (0.5 - diff(minMax(lat_lim))) / 2);
                nwse = [lat_lim(2), lon_lim(1), lat_lim(1), lon_lim(2)];
            else
                lon_lim = nwse([2 4]);
                lat_lim = nwse([3 1]);
                lon_lim = lon_lim + [-1 1] * max(0, (0.5 - diff(minMax(lon_lim))) / 2);
                lat_lim = lat_lim + [-1 1] * max(0, (0.5 - diff(minMax(lat_lim))) / 2);
            end
            clon = nwse([2 4]) + [-0.001 0.001];
            clat = max(-90, min(90, nwse([3 1]) + [-0.001 0.001]));

            if flag_dtm == 2
                subplot(1,2,1);
            end
            
            for i = 1 : iif(flag_dtm == 2, 2, 1)
                if flag_dtm == 2
                    subplot(1,2,i);
                end
                %m_proj('equidistant','lon',clon,'lat',clat);   % Projection
                m_proj('utm', 'lon',lon_lim,'lat',lat_lim);   % Projection
                axes
                cmap = flipud(gray(1000)); colormap(cmap(150: end, :));
                drawnow;
            end
            
            % retrieve external DTM
            try
                [dtm, lat_dtm, lon_dtm] = Core.getRefDTM(nwse, 'ortho', 'low');
                dtm = flipud(dtm);
                % comment the following line to have bathimetry
                dtm(dtm < 0) = nan; % - 1/3 * max(dtm(:));
                
                % uncomment the following line to have colors
                % colormap(Cmap.adaptiveTerrain(minMax(dtm(:))));
                % drawnow;
                
                [shaded_dtm, x, y] = m_shadedrelief(lon_dtm, lat_dtm, dtm, 'nan', [0.98, 0.98, 1], 'gra', 30);
                %h_dtm = m_pcolor(lon_dtm, lat_dtm, dtm);
                %h_dtm.CData = shaded_dtm;
                for i = 1 : iif(flag_dtm == 2, 2, 1)
                    if flag_dtm == 2
                        subplot(1,2,i);
                    end
                    m_image(lon_dtm, lat_dtm, shaded_dtm);
                end
            catch
                % use ETOPO1 instead
                for i = 1 : iif(flag_dtm == 2, 2, 1)
                    if flag_dtm == 2
                        subplot(1,2,i);
                    end
                    m_etopo2('shadedrelief','gradient', 3);
                end
            end
            
            for i = 1 : iif(flag_dtm == 2, 2, 1)
                if flag_dtm == 2
                    subplot(1,2,i);
                end
                % read shapefile
                shape = 'none';
                if (~strcmp(shape,'none'))
                    if (~strcmp(shape,'coast')) && (~strcmp(shape,'fill'))
                        if (strcmp(shape,'10m'))
                            M = m_shaperead('countries_10m');
                        elseif (strcmp(shape,'30m'))
                            M = m_shaperead('countries_30m');
                        else
                            M = m_shaperead('countries_50m');
                        end
                        [x_min, y_min] = m_ll2xy(min(lon_lim), min(lat_lim));
                        [x_max, y_max] = m_ll2xy(max(lon_lim), max(lat_lim));
                        for k = 1 : length(M.ncst)
                            lam_c = M.ncst{k}(:,1);
                            ids = lam_c <  min(lon);
                            lam_c(ids) = lam_c(ids) + 360;
                            phi_c = M.ncst{k}(:,2);
                            [x, y] = m_ll2xy(lam_c, phi_c);
                            if sum(~isnan(x))>1
                                x(find(abs(diff(x)) >= abs(x_max - x_min) * 0.90) + 1) = nan; % Remove lines that occupy more than th 90% of the plot
                                line(x,y,'color', [0.3 0.3 0.3]);
                            end
                        end
                    else
                        if (strcmp(shape,'coast'))
                            m_coast('line','color', lineCol);
                        else
                            m_coast('patch',lineCol);
                        end
                    end
                end            
                hold on;
                
                % Enable box
                m_grid('box','fancy','tickdir','in', 'fontsize', 16);
                % m_ruler(1.1, [.05 .40], 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
                drawnow
                m_ruler([.7 1], -0.08, 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
            end
            
            [tropo_grid, x_grid, y_grid, time, tropo_height_correction, tropo_clim] = sta_list.getTropoMap(par_name, rate);
            if flag_dtm == 1
                tropo_grid = tropo_grid + repmat(tropo_height_correction,1,1,size(tropo_grid,3));
            end
            
            if flag_dtm == 2
                ax = subplot(1,2,1);
            end
            imh = m_pcolor(x_grid, y_grid, tropo_grid(:,:,1));
            imh.FaceAlpha = 0.95;
            caxis(tropo_clim(1,:));

            switch lower(par_name)
                case 'ztd'
                case 'zwd'
                    caxis([max(0, tropo_clim(1,1)) min(48, tropo_clim(1,2))]);
                    %tropo_clim = caxis()]
                case 'gn'
                case 'ge'
                case 'pwv'
                case 'zhd'
                case 'nsat'
            end
            %cmap = Cmap.get('c51',512);
            %colormap(flipud(cmap(2:end,:)));
            cmap = Cmap.get('RdBu'); 
            colormap(cmap);
            %colormap(Cmap.smoothMap(Cmap.noaaRain));
            
            % redraw boxes
            m_grid('box','fancy','tickdir','in', 'fontsize', 16);
            %colormap(flipud(gat(1024, false)));
            %colormap(Cmap.get('viridis', 32));
            if flag_dtm == 2  
                cax = m_contfbar([.05 .55], 0.04, tropo_clim(1, 1), tropo_clim(1) : (diff(tropo_clim(1,:)) / size(cmap,1)) : tropo_clim(1, 2) ,'edgecolor','none','endpiece','no', 'fontsize', 16);
                xlabel(cax,'cm','color','k');
            else
                cax = m_contfbar([.15 .55], -0.05, tropo_clim(1, 1), tropo_clim(1) : (diff(tropo_clim(1,:)) / size(cmap,1)) : tropo_clim(1, 2) ,'edgecolor','none','endpiece','no', 'fontsize', 16);
                xlabel(cax,'cm','color','k');
            end
            xlabel(cax,'cm','color','k');
            
            % ax = fig_handle.Children(end);
            % ax.FontSize = 20;
            % ax.FontWeight = 'bold';
            % if new_fig
            %     if FTP_Downloader.checkNet()
            %         plot_google_map('alpha', 0.65, 'MapType', 'satellite');
            %     end
            %     xlabel('Longitude [deg]');
            %     ylabel('Latitude [deg]');
            % end
            
            if flag_dtm == 1
                th = title(sprintf([par_str ' map @%s at ground level\\fontsize{5} \n'], time.getEpoch(1).toString('yyyy-mm-dd HH:MM:SS')), 'FontSize', 20);
            else
                th = title(sprintf([par_str ' map @%s at sea level\\fontsize{5} \n'], time.getEpoch(1).toString('yyyy-mm-dd HH:MM:SS')), 'FontSize', 20);
            end
            
            if flag_dtm == 2  
                %tropo_clim(2,:) = [max(0, tropo_clim(1) + min(tropo_height_correction(:))) (tropo_clim(2) + max(tropo_height_correction(:)))];
                % uniform axes
                ax2 = subplot(1,2,2);
                imh2 = m_pcolor(x_grid, y_grid, tropo_grid(:,:,1) + tropo_height_correction);
                imh2.FaceAlpha = 0.95;
                caxis(tropo_clim(2,:)); 
                %colormap(Core_UI.CMAP_51(2:end,:));
                %colormap(flipud(gat(1024, false)));
                %colormap(gat2);       
                
                switch lower(par_name)
                    case 'ztd'
                    case 'zwd'
                        caxis([max(0, tropo_clim(1,1)) min(48, tropo_clim(1,2))]);
                        %tropo_clim = [tropo_clim(1,:); caxis()]; 
                    case 'gn'
                    case 'ge'
                    case 'pwv'
                    case 'zhd'
                    case 'nsat'
                end
                %cmap = Cmap.get('c51',512);
                cmap = Cmap.get('RdBu');
                colormap(flipud(cmap(2:end,:)));
                %colormap(Cmap.smoothMap(Cmap.noaaRain));
            
                % redraw boxes
                m_grid('box','fancy','tickdir','in', 'fontsize', 16);
                % ax2.FontSize = 20;
                % ax2.FontWeight = 'bold';
                % if new_fig
                %    if FTP_Downloader.checkNet()
                %        plot_google_map('alpha', 0.65, 'MapType', 'satellite');
                %    end
                %    xlabel('Longitude [deg]');
                %    ylabel('Latitude [deg]');
                %end
                if flag_dtm == 2
                    cax2 = m_contfbar([.05 .55], 0.04, tropo_clim(1, 1), tropo_clim(1) : (diff(tropo_clim(1,:)) / size(cmap,1)) : tropo_clim(1, 2) ,'edgecolor','none','endpiece','no', 'fontsize', 16);
                    xlabel(cax2,'cm','color','k');
                else
                    cax2 = m_contfbar([.15 .55], -0.05, tropo_clim(1, 1), tropo_clim(1) : (diff(tropo_clim(1,:)) / size(cmap,1)) : tropo_clim(1, 2) ,'edgecolor','none','endpiece','no', 'fontsize', 16);
                    xlabel(cax2,'cm','color','k');
                end
                th2 = title(ax2, sprintf('at ground level\\fontsize{5} \n'), 'FontSize', 20);
            end
            
            f.Visible = 'on'; drawnow;
            %Core.getLogger.addMarkedMessage('Press any key to start playing');
            %pause

            % Add logos
            if flag_export
                fprintf('Stopped in debug mode to check the export size and position of the elements\nType dbcont to continue ')
                % fh = gcf; fh.Units = 'pixel'; fh.Position([3 4]) = [1100 960];
                keyboard % to check before export that everything is aligned
                if flag_dtm ~= 2
                    fh = gcf; fh.Units = 'pixel'; fh.Position([3 4]) = [1040 764];
                    drawnow
                    drawnow
                end
                Core_UI.insertLogo(f, 'SouthEast');
                warning off;                
                im = {};
                f.Visible = 'off';
                Core.getLogger.addMarkedMessage('Exporting video');
                fprintf('%5d/%5d', 0, 99999);
                
                file_name = sprintf('AniMap%sInterp-%s-%s-at%s', par_str_short, time.first.toString('yyyymmdd_HHMMSS'), time.last.toString('yyyymmdd_HHMMSS'), datestr(now, 'yyyymmdd_HHMMSS'));
                if ismac() || ispc()
                    % Better compression on Mac > 10.7 and Win > 7
                    video_out = VideoWriter(fullfile(Core.getState.getOutDir, [file_name '.mp4']), 'MPEG-4');
                else
                    % Linux doesn't have mp4 compression avaiable
                    video_out = VideoWriter(fullfile(Core.getState.getOutDir, [file_name '.avi']));
                end
                video_out.FrameRate = 10;
                video_out.Quality = 91;
                open(video_out);
            else
                Core_UI.insertLogo(f, 'SouthEast');
                f.Visible = 'on'; drawnow;
            end
            drawnow
            
            for i = 1 : time.length
                if any(serialize(tropo_grid(:,:,i)))
                    th.String = sprintf([par_str ' map %s at sea level'], time.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS'));
                    imh.CData = tropo_grid(:,:,i);
                    %imh.AlphaData = ~isnan(tropo_grid(:,:,i));
                    if flag_dtm == 2
                        imh2.CData = tropo_grid(:,:,i) + tropo_height_correction;
                        %imh2.AlphaData = imh.AlphaData;
                    end
                    if not(flag_export)
                        drawnow
                    end
                end
                if flag_export
                    fprintf('%s%5d/%5d',char(8 * ones(1,11)), i, time.length);
                    frame = getframe(f);
                    ss = 1; % subsample (1:2)
                    writeVideo(video_out, frame(1:ss:end,1:ss:end,:)); 
                end
            end
            
            if flag_export
                fprintf('%s',char(8 * ones(1,11)));
                close(video_out);
                Core.getLogger.addStatusOk(sprintf('"%s" done ^_^', fullfile(Core.getState.getOutDir, video_out.Filename)));
                close(f);
                warning on;
            end
        end

        function fh_list = showAniMapTropoInterpAndDiff(sta_list_full, par_name, nwse, rate, flag_dtm, flag_export)
            % Show a tropo map with all the station in sta_list
            %
            % INPUT
            %   tropo_par   accepted tropo parameter:
            %               - 'zwd'
            %               - 'ztd'
            %   nwse        Nort West South East coordinates
            %   rate        rate in seconds, nearest to closest observation 
            %               it should be a subsample of the data rate (e.g. 300 with 30s data)
            %   flag_dtm    flag to add height_correction (default == false)
            %   flag_export if true try to export a video (the frames are not going to be seen)
            %               the video file will be saved in the out folder specified in the project
            %
            % SYNTAX
            %   sta_list.showAniMapTropoInterpAndDiff(par_name, <new_fig>, <rate>, <flag_dtm>, <flag_export>);
            %
            % EXAMPLE
            %   % over Japan
            %   sta_list.showAniMapTropoInterpAndDiff('ZWD', [45.8, 123.5, 23, 146.5], 200, 2, false);
            
            switch lower(par_name)
                case 'ztd'
                    par_str = 'ZTD';
                    par_str_short = 'ZTD';
                case 'zwd'
                    par_str = 'ZWD';
                    par_str_short = 'ZWD';
                case 'gn'
                case 'ge'
                case 'pwv'
                    par_str = 'PWV';
                    par_str_short = 'PWV';
                case 'zhd'
                case 'nsat'
            end                        

            f = figure;
            
            fh_list = f;
            fig_name = sprintf('RecAniMapInterpDiff_%s', upper(par_name));
            f.UserData = struct('fig_name', fig_name);
            
            if nargin < 4 || isempty(rate)
                rate = 300;
            end
            if nargin < 5
                flag_dtm = false;
            end
            if nargin < 6
                flag_export = false;
            end

            maximizeFig(f);
            
            f.Visible = 'off';
            f.Color = [1 1 1];

            % Set map projection / limits
            
            [lat, lon] = sta_list_full.getMedianPosGeodetic();
            if nargin >= 3 && ~isempty(nwse)
                margin = 0.5;
                id_ok = (lat >= (nwse(3) - margin)) & (lat <= (nwse(1) + margin)) & ...
                    (lon >= (nwse(2) - margin)) & (lon <= (nwse(4) + margin));
                sta_list = sta_list_full(id_ok);
            end
            
            if nargin < 3 || isempty(nwse)
                sta_list = sta_list_full;
                % set map limits
                if numel(sta_list) == 1
                    lon_lim = minMax(lon) + [-0.05 0.05];
                    lat_lim = minMax(lat) + [-0.05 0.05];
                else
                    lon_lim = minMax(lon); lon_lim = lon_lim + [-1 1] * diff(lon_lim)/15;
                    lat_lim = minMax(lat); lat_lim = lat_lim + [-1 1] * diff(lat_lim)/15;
                end
                lon_lim = lon_lim + [-1 1] * max(0, (0.5 - diff(minMax(lon_lim))) / 2);
                lat_lim = lat_lim + [-1 1] * max(0, (0.5 - diff(minMax(lat_lim))) / 2);
                nwse = [lat_lim(2), lon_lim(1), lat_lim(1), lon_lim(2)];
            else
                lon_lim = nwse([2 4]);
                lat_lim = nwse([3 1]);
                lon_lim = lon_lim + [-1 1] * max(0, (0.5 - diff(minMax(lon_lim))) / 2);
                lat_lim = lat_lim + [-1 1] * max(0, (0.5 - diff(minMax(lat_lim))) / 2);
            end
            clon = nwse([2 4]) + [-0.001 0.001];
            clat = max(-90, min(90, nwse([3 1]) + [-0.001 0.001]));

            if flag_dtm == 2
                subplot(1,2,1);
            end
            
            for i = 1 : iif(flag_dtm == 2, 2, 1)
                if flag_dtm == 2
                    subplot(1,2,i);
                end
                %m_proj('equidistant','lon',clon,'lat',clat);   % Projection
                m_proj('utm', 'lon',lon_lim,'lat',lat_lim);   % Projection
                axes
                cmap = flipud(gray(1000)); colormap(cmap(150: end, :));
                drawnow;
            end
            
            % retrieve external DTM
            try
                [dtm, lat_dtm, lon_dtm] = Core.getRefDTM(nwse, 'ortho', 'low');
                dtm = flipud(dtm);
                % comment the following line to have bathimetry
                dtm(dtm < 0) = nan; % - 1/3 * max(dtm(:));
                
                % uncomment the following line to have colors
                % colormap(Cmap.adaptiveTerrain(minMax(dtm(:))));
                % drawnow;
                
                [shaded_dtm, x, y] = m_shadedrelief(lon_dtm, lat_dtm, dtm, 'nan', [0.98, 0.98, 1], 'gra', 30);
                %h_dtm = m_pcolor(lon_dtm, lat_dtm, dtm);
                %h_dtm.CData = shaded_dtm;
                for i = 1 : iif(flag_dtm == 2, 2, 1)
                    if flag_dtm == 2
                        subplot(1,2,i);
                    end
                    m_image(lon_dtm, lat_dtm, shaded_dtm);
                end
            catch
                % use ETOPO1 instead
                for i = 1 : iif(flag_dtm == 2, 2, 1)
                    if flag_dtm == 2
                        subplot(1,2,i);
                    end
                    m_etopo2('shadedrelief','gradient', 3);
                end
            end
            
            for i = 1 : iif(flag_dtm == 2, 2, 1)
                if flag_dtm == 2
                    subplot(1,2,i);
                end
                % read shapefile
                shape = 'none';
                if (~strcmp(shape,'none'))
                    if (~strcmp(shape,'coast')) && (~strcmp(shape,'fill'))
                        if (strcmp(shape,'10m'))
                            M = m_shaperead('countries_10m');
                        elseif (strcmp(shape,'30m'))
                            M = m_shaperead('countries_30m');
                        else
                            M = m_shaperead('countries_50m');
                        end
                        [x_min, y_min] = m_ll2xy(min(lon_lim), min(lat_lim));
                        [x_max, y_max] = m_ll2xy(max(lon_lim), max(lat_lim));
                        for k = 1 : length(M.ncst)
                            lam_c = M.ncst{k}(:,1);
                            ids = lam_c <  min(lon);
                            lam_c(ids) = lam_c(ids) + 360;
                            phi_c = M.ncst{k}(:,2);
                            [x, y] = m_ll2xy(lam_c, phi_c);
                            if sum(~isnan(x))>1
                                x(find(abs(diff(x)) >= abs(x_max - x_min) * 0.90) + 1) = nan; % Remove lines that occupy more than th 90% of the plot
                                line(x,y,'color', [0.3 0.3 0.3]);
                            end
                        end
                    else
                        if (strcmp(shape,'coast'))
                            m_coast('line','color', lineCol);
                        else
                            m_coast('patch',lineCol);
                        end
                    end
                end            
                hold on;
                
                % Enable box
                m_grid('box','fancy','tickdir','in', 'fontsize', 16);
                % m_ruler(1.1, [.05 .40], 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
                drawnow
                m_ruler([.7 1], -0.08, 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
            end
                        
            [tropo_grid, x_grid, y_grid, time, tropo_height_correction, tropo_clim] = sta_list.getTropoMap(par_name, rate);
            if flag_dtm == 1
                tropo_grid = tropo_grid + tropo_height_correction;
            end
            
            if flag_dtm == 2
               ax = subplot(1,2,1);
            else
                ax = gca;
            end
            imh = m_pcolor(x_grid, y_grid, tropo_grid(:,:,1));
            imh.FaceAlpha = 0.95;
            caxis(tropo_clim(1,:));
            
            switch lower(par_name)
                case 'ztd'
                case 'zwd'
                    caxis([max(0, tropo_clim(1,1)) min(48, tropo_clim(1,2))]);
                    %tropo_clim = caxis()];
                case 'gn'
                case 'ge'
                case 'pwv'
                case 'zhd'
                case 'nsat'
            end
            cmap = Cmap.get('c51',512);
            colormap(flipud(cmap(2:end,:)));
            %colormap(Cmap.smoothMap(Cmap.noaaRain));

            % redraw boxes
            m_grid('box','fancy','tickdir','in', 'fontsize', 16);
            %colormap(flipud(gat(1024, false)));
            %colormap(Cmap.get('viridis', 32));
            th = title(sprintf([par_str ' map @%s at sea level\\fontsize{5} \n'], time.getEpoch(1).toString('yyyy-mm-dd HH:MM:SS')), 'FontSize', 20);
            if flag_dtm > 0
                tropo_diff = (tropo_grid(:,:,2 : end) - tropo_grid(:,:,1 : end - 1)) / time.getRate * 3600;
                tropo_diff_clim = [-1 1] * perc(abs(serialize(tropo_diff)), 0.998);
                % remove outliers
                tropo_diff(tropo_diff > 2 * max(tropo_diff_clim)) = nan;
                tropo_diff(tropo_diff < 2 * min(tropo_diff_clim)) = nan;
            end
            if flag_dtm == 2  
                cax = m_contfbar([.05 .55], -0.04, tropo_clim(1, 1), tropo_clim(1) : (diff(tropo_clim(1,:)) / size(cmap,1)) : tropo_clim(1, 2) ,'edgecolor','none','endpiece','no', 'fontsize', 16);
                xlabel(cax,'cm','color','k');                
            else
                if flag_dtm == 1
                    caxis(tropo_diff_clim);
                    colormap(ax,flipud(Cmap.get('RdBu')));
                    cax = m_contfbar([.05 .55], 0.04, tropo_diff_clim(1, 1), tropo_diff_clim(1) : (diff(tropo_diff_clim(1,:)) / size(cmap,1)) : tropo_diff_clim(1, 2) ,'edgecolor','none','endpiece','no', 'fontsize', 16);
                    xlabel(cax, 'cm / h','color','k');
                    th = title(ax, sprintf('Temporal derivative of ZWD\n@%s\\fontsize{5} \n', time.getEpoch(1).toString('yyyy-mm-dd HH:MM:SS')), 'FontSize', 20);
                else
                    cax = m_contfbar([.15 .55], -0.05, tropo_clim(1, 1), tropo_clim(1) : (diff(tropo_clim(1,:)) / size(cmap,1)) : tropo_clim(1, 2) ,'edgecolor','none','endpiece','no', 'fontsize', 16);
                    xlabel(cax,'cm','color','k');
                end
            end
            
                        
            if flag_dtm == 2  
                %tropo_clim(2,:) = [max(0, tropo_clim(1) + min(tropo_height_correction(:))) (tropo_clim(2) + max(tropo_height_correction(:)))];
                % uniform axes
                ax2 = subplot(1,2,2);
                imh2 = m_pcolor(x_grid, y_grid, nan(size(tropo_grid(:,:,1))));
                imh2.FaceAlpha = 1;
                %caxis(tropo_clim(2,:)); 
                %colormap(Core_UI.CMAP_51(2:end,:));
                %colormap(flipud(gat(1024, false)));
                %colormap(gat2);                            
                % redraw boxes
                m_grid('box','fancy','tickdir','in', 'fontsize', 16);
                % ax2.FontSize = 20;
                % ax2.FontWeight = 'bold';
                % if new_fig
                %    if FTP_Downloader.checkNet()
                %        plot_google_map('alpha', 0.65, 'MapType', 'satellite');
                %    end
                %    xlabel('Longitude [deg]');
                %    ylabel('Latitude [deg]');
                %end
                %tropo_diff_clim = [-1 1] * max(abs(minMax(serialize(tropo_grid(:,:,2 : end) - tropo_grid(:,:,1 : end - 1)))));
                if flag_dtm == 2
                    colormap(ax2,flipud(Cmap.get('RdBu')));
                    ax2.CLim = tropo_diff_clim;
                    cax2 = m_contfbar([.05 .55], 0.04, tropo_diff_clim(1, 1), tropo_diff_clim(1) : (diff(tropo_diff_clim(1,:)) / size(cmap,1)) : tropo_diff_clim(1, 2) ,'edgecolor','none','endpiece','no', 'fontsize', 16);
                    xlabel(cax2,'cm / h','color','k');
                    th2 = title(ax2, sprintf('Temporal derivative of ZWD\\fontsize{5} \n'), 'FontSize', 20);
                end
            end
            
            f.Visible = 'on'; drawnow;
            %Core.getLogger.addMarkedMessage('Press any key to start playing');
            %pause

            % Add logos
            if flag_export
                fprintf('Stopped in debug mode to check the export size and position of the elements\nType dbcont to continue ')
                % fh = gcf; fh.Units = 'pixel'; fh.Position([3 4]) = [1100 960];
                keyboard % to check before export that everything is aligned
                if flag_dtm ~= 2
                    %fh = gcf; fh.Units = 'pixel'; fh.Position([3 4]) = [1040 764];
                    drawnow
                    drawnow
                end
                Core_UI.insertLogo(f, 'NorthWest');
                warning off;                
                im = {};
                f.Visible = 'on'; drawnow;
                Core.getLogger.addMarkedMessage('Exporting video');
                fprintf('%5d/%5d', 0, 99999);
                
                file_name = sprintf('AniMap%sInterpAndDiff-%s-%s-at%s', par_str_short, time.first.toString('yyyymmdd_HHMMSS'), time.last.toString('yyyymmdd_HHMMSS'), datestr(now, 'yyyymmdd_HHMMSS'));
                if ismac() || ispc()
                    % Better compression on Mac > 10.7 and Win > 7
                    video_out = VideoWriter(fullfile(Core.getState.getOutDir, [file_name '.mp4']), 'MPEG-4');
                else
                    % Linux doesn't have mp4 compression avaiable
                    video_out = VideoWriter(fullfile(Core.getState.getOutDir, [file_name '.avi']));
                end
                
                video_out.FrameRate = 30;
                video_out.Quality = 91;
                open(video_out);
            else
                Core_UI.insertLogo(f, 'SouthEast');
                f.Visible = 'on'; drawnow;
            end
            drawnow
            
            % compute cut_mask (trim white borders)
            frame = getframe(f);
            id_x = find(any(sum(frame.cdata,3) ~= 765), 1, 'first'):find(any(sum(frame.cdata,3) ~= 765), 1, 'last');
            id_y = find(any(sum(frame.cdata,3) ~= 765, 2), 1, 'first'):find(any(sum(frame.cdata,3) ~= 765, 2), 1, 'last');
            % end of trim
            
            for i = 2 : time.length
                if any(serialize(tropo_grid(:,:,i)))
                    if flag_dtm == 1
                        th = title(ax, sprintf('Temporal derivative of ZWD at sea level\n@%s\\fontsize{5} \n', time.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS')), 'FontSize', 20);
                    else
                        th.String = sprintf([par_str ' map %s at sea level'], time.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS'));
                    end
                    if flag_dtm == 1
                        imh.CData = Cmap.apply(tropo_diff(:,:,i-1), flipud(Cmap.get('RdBu')), tropo_diff_clim, [nan, nan, nan]);
                    else
                        imh.CData = tropo_grid(:,:,i);
                    end
                    %imh.AlphaData = ~isnan(tropo_grid(:,:,i));
                    if flag_dtm == 2
                        imh2.CData = Cmap.apply(tropo_diff(:,:,i-1), flipud(Cmap.get('RdBu')), tropo_diff_clim, [nan, nan, nan]);
                        %imh2.FaceAlpha = 1;
                        %imh2.CData = tropo_grid(:,:,i) - tropo_grid(:,:,i-1);
                        %imh2.AlphaData = imh.AlphaData;
                    end
                    if not(flag_export)
                        drawnow
                    end
                end
                if flag_export
                    fprintf('%s%5d/%5d',char(8 * ones(1,11)), i, time.length);
                    frame = getframe(f);
                    frame.cdata = frame.cdata(id_y,id_x,:);
                    writeVideo(video_out, frame); 
                end
            end
            
            if flag_export
                fprintf('%s',char(8 * ones(1,11)));
                close(video_out);
                Core.getLogger.addStatusOk(sprintf('"%s" done ^_^', fullfile(Core.getState.getOutDir, video_out.Filename)));
                close(f);
                warning on;
            end
        end

        function fh_list = showDt(this)
            % Plot Clock error
            %
            % SYNTAX
            %   sta_list.plotDt

            fh_list = [];
            for r = 1 : size(this, 2)
                rec = this(r);
                if ~isempty(rec)
                    fh_list = [fh_list; rec.out.showDt()]; %#ok<AGROW>
                end
            end
        end
        
        function fh_list = showOutliersAndCycleSlip(sta_list, sys_list)
            % Show Outliers and cycle slips for each receiver workspace
            % (cartesian plot)
            %
            % SYNTAX
            %   this.showOutliersAndCycleSlip(sys_list)
            
            fh_list = [];
            for s = 1 : numel(sta_list)
                if nargin < 2
                    fh_list = [fh_list; sta_list(s).work.showOutliersAndCycleSlip()]; %#ok<AGROW>
                else
                    fh_list = [fh_list; sta_list(s).work.showOutliersAndCycleSlip(sys_list)]; %#ok<AGROW>
                end
            end
        end
        
        function fh_list = showOutliersAndCycleSlip_p(sta_list, sys_list)
            % Show Outliers and cycle slips for each receiver workspace
            % (polar plot)
            %
            % SYNTAX
            %   this.showOutliersAndCycleSlip(sys_list)
            
            if nargin < 2
                sys_list = Core.getConstellationCollector.getAvailableSys();
            end
            fh_list = [];
            for s = 1 : numel(sta_list)
                fh_list = [fh_list; sta_list(s).work.showOutliersAndCycleSlip_p(sys_list)]; %#ok<AGROW>
            end
        end
        
        function fh_list = showResSkyCartScatter(sta_list, sys_list, type)
            % Show Residuals for each receiver workspace
            % (cartesian plot)
            %
            % SYNTAX
            %   this.showResSkyCartScatter(sys_list)
            
            if nargin < 2 || isempty()
                sys_list = Core.getConstellationCollector.getAvailableSys();
            end
            
            fh_list = [];
            for s = 1 : numel(sta_list)
                fh_list = [fh_list; sta_list(s).work.showResSkyCartScatter(sys_list)]; %#ok<AGROW>
            end
        end
        
        function fh_list = showResSkyPolarScatter(sta_list, sys_list)
            % Show Residuals for each receiver workspace
            % (polar plot)
            %
            % SYNTAX
            %   this.showResSkyPolarScatter(sys_list)
            
            if nargin < 2
                sys_list = Core.getConstellationCollector.getAvailableSys();
            end
            fh_list = [];
            for s = 1 : numel(sta_list)
                fh_list = [fh_list; sta_list(s).work.showResSkyPolarScatter(sys_list)]; %#ok<AGROW>
            end
        end
        
        function fh_list = showSNR_p(sta_list, sys_list)
            % Show SNR for each receiver workspace
            % (polar plot)
            %
            % SYNTAX
            %   this.showSNR_p(sys_list)
            
            if nargin < 2
                sys_list = Core.getConstellationCollector.getAvailableSys();
            end
            fh_list = [];
            for s = 1 : numel(sta_list)
                fh_list = [fh_list; sta_list(s).work.showSNR_p(sys_list)]; %#ok<AGROW>
            end
        end
        
        function fh_list = showSNR_z(sta_list, sys_list, l_max)
            % Show SNR for each receiver workspace
            % (polar plot Zernike interpolated)
            %
            % SYNTAX
            %   this.showSNR_p(sys_list)
            
            if nargin < 2 || isempty(sys_list)
                sys_list = Core.getConstellationCollector.getAvailableSys();
            end
            fh_list = [];
            for s = 1 : numel(sta_list)
                if nargin == 3
                    fh_list = [fh_list; sta_list(s).work.showSNR_z(sys_list, l_max)]; %#ok<AGROW>
                else
                    fh_list = [fh_list; sta_list(s).work.showSNR_z(sys_list)]; %#ok<AGROW>
                end
            end
        end
        
        function fh_list = showMultiPathModel(sta_list, mp_type, ant_mp)
            % Show MultiPath Maps for each receiver workspace
            % (polar plot Zernike interpolated)
            %
            % SYNTAX
            %   this.showSNR_p(sys_list)
            if nargin == 1
                mp_type = 1;
            end
            log = Core.getLogger;
            state = Core.getState;
            flag_show = 0;
            fh_list = [];
            for r = 1 : numel(sta_list)
                rec = sta_list(r);
                if nargin < 3 || isempty(ant_mp)
                    ant_mp = rec.getAntennaMultiPath();
                end
                if ~isempty(ant_mp)
                    % Get the satellite systems available in the zerniche multipath struct
                    sys_c_list = cell2mat(fields(ant_mp)');
                    for sys_c = sys_c_list
                        if isfield(ant_mp, sys_c)
                            % This constellation is already present into the applied Zernike MultiPath set of coefficients
                            trk_list = fields(ant_mp.(sys_c))';
                            for trk = trk_list
                                trk = strrep(trk{1}, ' ', '_'); %#ok<FXSET> % structures do not support spaces
                                if isfield(ant_mp.(sys_c), trk)
                                    switch (mp_type)
                                        case 1 % Zernike map
                                            str_type = 'Zernike ';
                                            str_out_type = 'Z_';
                                            mp_map = double(ant_mp.(sys_c).(trk).z_map);
                                        case 2 % Zernike map + gridded residuals
                                            str_type = 'Zernike + congruent ';
                                            str_out_type = 'ZC_';
                                            mp_map = double(ant_mp.(sys_c).(trk).r_map);
                                        case 3 % Simple Gridding of size [stk_grid_step]
                                            str_type = 'regular smoothed ';
                                            str_out_type = 'R_';
                                            mp_map = double(ant_mp.(sys_c).(trk).g_map);
                                        case 4 % Congruent cells gridding of size [stk_grid_step]
                                            str_type = 'congruent smoothed ';
                                            str_out_type = 'C_';
                                            mp_map = double(ant_mp.(sys_c).(trk).c_map);
                                        case 5 % Simple Gridding of size [1x1]
                                            str_type = 'regular ';
                                            str_out_type = 'G1_';
                                            mp_map = double(ant_mp.(sys_c).(trk).g1_map);
                                        case 6 % c1_map Congruent cells gridding of size [1x1]
                                            str_type = 'congruent ';
                                            str_out_type = 'C1_';
                                            mp_map = double(ant_mp.(sys_c).(trk).c1_map);
                                    end
                                    if numel(mp_map) <= 1
                                        log.addWarning(sprintf('The %c%s %s map is not available', sys_c, trk, str_type));
                                    else
                                        [az_grid, el_grid] = Core_Utils.getPolarGrid(360 / size(mp_map, 2), 90 / size(mp_map, 1));
                                        az_grid = Core_Utils.deg2rad(az_grid);
                                        el_grid = Core_Utils.deg2rad(el_grid);
                                        fh = figure('Visible', 'off'); fh.Name = sprintf('%03d: MP %s %s%s %s', fh.Number, rec.getMarkerName4Ch, sys_c, trk, str_type); fh.NumberTitle = 'off';
                                        polarImagesc(az_grid, (pi/2 - el_grid), 1e3*(mp_map));
                                        fh_list = [fh_list; fh];
                                        caxis(max(abs(minMax(caxis))) * [-1 1]);
                                        colormap((Cmap.get('PuOr', 2^11)));
                                        colorbar;
                                        fig_name = sprintf('MP_%sMap_%s_%s%s_%s', str_out_type, rec.getMarkerName4Ch, sys_c, trk, rec.getTime.last.toString('yyyymmdd_HHMM'));
                                        fh.UserData = struct('fig_name', fig_name);
                                        Core_UI.beautifyFig(fh, 'light');
                                        Core_UI.addExportMenu(fh);
                                        Core_UI.addBeautifyMenu(fh);
                                        fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                                        title((sprintf('Multipath %smitigation map of %s %s%s [mm]', str_type, rec.getMarkerName4Ch, sys_c, trk)), 'interpreter', 'none'); drawnow
                                        % caxis([-15 15]);
                                    end
                                end
                            end
                        end
                    end   
                    flag_show = flag_show + 1;
                end
            end
            if flag_show == 0 && numel(sta_list) > 0
                log.addError('Multipath mitigation model maps, no models found!');
            else
                log.addStatusOk(sprintf('%d of %d multipath mitigation models exported successfully', flag_show, numel(sta_list)));
            end
        end
        
        function fh_list = showQuality_p(sta_list, type, flag_smooth)
            % Plot Signal to Noise Ration in a skyplot
            % SYNTAX f_handles = this.plotSNR(sys_c)

            % SNRs
            if nargin < 2
                type = 'snr';
            end

            fh_list = [];

            for r = 1 : numel(sta_list)
                if ~sta_list(r).out.isEmpty
                    rec = sta_list(r).out;
                    [quality, az, el] = rec.getQuality();
                else
                    rec = sta_list(r).work;
                    [quality, az, el] = rec.getQuality(type);
                end

                if nargin > 2 && flag_smooth
                    quality = Receiver_Commons.smoothSatData([],[],zero2nan(quality), [], 'spline', 900 / this.getRate, 10); % smoothing Quality => to be improved
                end

                if (numel(az) ~= numel(quality))
                    log = Core.getLogger();
                    log.addError('Number of elements for az different from quality data\nPlotting id not possible');
                else
                    fh = figure('Visible', 'off'); fh.Name = sprintf('%03d: %s %s', fh.Number, upper(type), sta_list(r).getMarkerName4Ch); fh.NumberTitle = 'off';
                    
                    fh_list = [fh_list; fh]; %#ok<AGROW>
                    fig_name = sprintf('Quality_Polar_%s_%s', sta_list(r).getMarkerName4Ch, sta_list(r).getTime.first.toString('yyyymmdd_HHMM'));
                    fh.UserData = struct('fig_name', fig_name);

                    id_ok = (~isnan(quality));
                    polarScatter(serialize(az(id_ok)) / 180 * pi, serialize(90 - el(id_ok)) / 180 * pi, 45, serialize(quality(id_ok)), 'filled');
                    colormap(jet);  cax = caxis();
                    switch type
                        case 'snr'
                            caxis([min(cax(1), 10), max(cax(2), 55)]);
                            setColorMap('jet', [10 55], 0.9);
                    end
                    colorbar();
                    h = title(sprintf('%s - receiver %s', upper(type), sta_list(r).getMarkerName4Ch()), 'interpreter', 'none');
                    h.FontWeight = 'bold'; h.Units = 'pixels';
                    
                    Core_UI.beautifyFig(fh);
                    Core_UI.addExportMenu(fh);
                    Core_UI.addBeautifyMenu(fh);
                    fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                end
            end
            
        end

        function fh_list = showResPerSat(sta_list, sys_c_list, type)
                        % Plot the residuals of phase per Satellite
            %
            % INPUT
            %   type    can be:
            %            'co'   -> Combined residuals (one set for each satellite) DEFAULT
            %            'pr'   -> Uncombined pseudo-ranges residuals
            %            'ph'   -> Uncombined carrier-phase residuals
            %   res     is the matrix of residuals satellite by satellite and can be passed from e.g. NET
            %
            % SYNTAX
            %   this.showResPerSat(sys_c_list, res, type)

            if nargin < 2 || isempty(sys_c_list)
                sys_c_list = Core.getConstellationCollector.getAvailableSys();
            end
            if nargin < 3 || isempty(type)
                type = 'co';
            end
            
            fh_list = [];
            for r = 1 : size(sta_list, 2)
                rec = sta_list(r);
                if ~isempty(rec)
                    if ~rec.out.isEmpty
                        fh_list = [fh_list; rec.out.showResPerSat(sys_c_list, type)]; %#ok<AGROW>
                    else
                        fh_list = [fh_list; rec.work.showResPerSat(sys_c_list, type)]; %#ok<AGROW>
                    end
                end
            end
        end

        function fh_list = showRes(sta_list)
            % Plot Satellite Residuals
            %
            % SYNTAX
            %   sta_list.showRes()

            fh_list = [];
            for r = 1 : size(sta_list, 2)
                rec = sta_list(r);
                if ~isempty(rec)
                    if ~rec.out.isEmpty
                        fh_list = [fh_list; rec.out.showRes()]; %#ok<AGROW>
                    else
                        fh_list = [fh_list; rec.work.showRes()]; %#ok<AGROW>
                    end
                end
            end
        end

        function fh_list = showResMap(sta_list, step, sys_c_list, mode)
            % Plot Satellite Residuals as a map
            %
            % SYNTAX
            %   sta_list.showResMap(step)
            fh_list = [];
            if nargin < 2 || isempty(step)
                step = 0.5;
            end
            for r = 1 : size(sta_list, 2)
                rec = sta_list(r);
                if ~isempty(rec)
                    if ~rec.out.isEmpty
                        cc = rec.out.getCC;
                    else
                        cc = rec.work.getCC;
                    end
                    if nargin < 3 || isempty(sys_c_list)
                        sys_c_list = cc.getActiveSysChar;
                    end
                    for ss = sys_c_list
                        if ~rec.out.isEmpty
                            [map, map_fill, ~, az_g, el_g] = rec.out.getResMap(step, 3, ss);
                        else
                            [map, map_fill, ~, az_g, el_g] = rec.work.getResMap(step, 3, ss);
                        end
                        % restore the original mean data where observations are present
                        %map_fill(~isnan(map)) = map(~isnan(map));
                        
                        if ~any(map(:))
                            Core.getLogger.addWarning(sprintf('No data found for %s@%c', rec.getMarkerName4Ch, ss));
                        else
                            
                            f = figure('Visible', 'off');
                            
                            fh_list = [fh_list; f]; %#ok<AGROW>
                            fig_name = sprintf('Res_map_%s_%s_%s', rec.getMarkerName4Ch, cc.getSysName(ss), rec.getTime.first.toString('yyyymmdd_HHMM'));
                            f.UserData = struct('fig_name', fig_name);
                            
                            f.Name = sprintf('%03d: ResMap %s@%c', f.Number, rec.getMarkerName4Ch, ss); f.NumberTitle = 'off';
                            if (nargin < 4) || isempty(mode)
                                mode = 'cart';
                            end
                            switch mode
                                case 'cart'
                                    %% Cartesian projection
                                    %img = imagesc(az_g, el_g, 1e3 * circshift(abs(map_fill), size(map_fill, 2) / 2, 2));
                                    img = imagesc(az_g, el_g, 1e3 * map_fill);
                                    set(gca,'YDir','normal');
                                    grid on
                                    image_alpha = 0.5; % everywhere 1 where obs are present
                                    %img.AlphaData = (~isnan(circshift(abs(map), size(map, 2) / 2, 2)) * 0.7) + 0.3;
                                    img.AlphaData = (~isnan(map) * (1 - image_alpha)) + image_alpha;
                                    %colormap(flipud(hot)); colorbar(); caxis([0, 0.02]);
                                    
                                    %caxis(1e3 * [min(abs(map(:))) min(20, min(6*std(zero2nan(map(:)),'omitnan'), max(abs(zero2nan(map(:))))))]);
                                    caxis(1e3 * perc(abs(map(:)), 0.99) * [-1 1]);
                                    colormap(Cmap.get('PuOr', 256));
                                    f.Color = [.95 .95 .95]; colorbar(); ax = gca; ax.Color = 'none';
                                    h = title(sprintf('Satellites residuals [mm] - receiver %s - %c', ss, rec.getMarkerName4Ch, ss),'interpreter', 'none');
                                    h.FontWeight = 'bold';
                                    hl = xlabel('Azimuth [deg]'); hl.FontWeight = 'bold';
                                    hl = ylabel('Elevation [deg]'); hl.FontWeight = 'bold';
                                    %ax = gca;
                                    %ax.PlotBoxAspectRatio = [1 1 1];
                                    %ax.DataAspectRatio(2) = ax.DataAspectRatio(1);
                                case '3D'
                                    %% 3D projection
                                    %                                 clf
                                    %                                 polarplot3d(1e3 * map_fill(1:2:end, 1:2:end),  'RadialRange',[-180 180] / 180 * pi, ...
                                    %                                     'AxisLocation', 0, 'InterpMethod', 'cubic', ...
                                    %                                     'PlotType', 'surfn', 'tickspacing', 15, ...
                                    %                                     'GridColor', [0.7 0.7 0.7]);
                                    polarplot3d(1e3 * map_fill, ...
                                        'AxisLocation', 0, 'InterpMethod', 'cubic', ...
                                        'PlotType', 'surfn', 'tickspacing', 15, ...
                                        'GridColor', [0.7 0.7 0.7]);
                                    caxis(1e3 * perc(abs(map(:)), 0.95) * [-1 1]);
                                    colormap(Cmap.get('RdGy', 256));
                                    colorbar();
                                    ax = gca;
                                    ax.PlotBoxAspectRatio = [1 1 1];
                                    ax.DataAspectRatio(2) = ax.DataAspectRatio(1);
                                    smap = ax.Children(end);
                                    smap.AlphaData = double(~isnan(map));
                                    smap.AlphaData(smap.AlphaData == 0) = smap.AlphaData(smap.AlphaData == 0) + 0.5;
                                    smap.FaceAlpha = 'flat';
                            end
                            Core_UI.beautifyFig(gcf);
                            Core_UI.addExportMenu(gcf);             
                            Core_UI.addBeautifyMenu(gcf);             
                            f.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                        end
                    end
                end
            end
        end

        function fh_list = showZtdSlant(sta_list, time_start, time_stop)
            fh_list = [];
            for r = 1 : size(sta_list, 2)
                rec = sta_list(~sta_list(r).isEmpty, r);
                if isempty(rec)
                    log = Core.getLogger();
                    log.addWarning('ZTD and/or slants have not been computed');
                else
                    if nargin < 3
                        fh_list = [fh_list; rec.out.showZtdSlant()]; %#ok<AGROW>
                    else
                        fh_list = [fh_list; rec.out.showZtdSlant(time_start, time_stop)]; %#ok<AGROW>
                    end
                end
            end
        end

        function fh_list = showPTH(sta_list)
            % Show plots for pressure, temperature and humidity
            %
            % SYNATAX
            %   sta_list.showPTH()
            [pressure, temperature, humidity, p_time, id_sync] = sta_list.getPTH_mr();

            f = figure('Visible', 'off');
            
            fh_list = f;
            if numel(sta_list) == 1
                % If I have only one receiver use as name the name of the receiver
                fig_name = sprintf('%s_%s_%s', 'PTH', sta_list.getMarkerName4Ch, sta_list.getTime.first.toString('yyyymmdd_HHMM'));
            else
                % If I have more than one receiver use as name the name of the project
                fig_name = sprintf('%s_%s', 'PTH', strrep(Core.getState.getPrjName,' ', '_'));
            end
            f.UserData = struct('fig_name', fig_name);
            
            f.Name = sprintf('%03d: %s %s', f.Number, 'PTH', sta_list(1).out.getCC.sys_c); f.NumberTitle = 'off';
            set(f,'defaultAxesColorOrder', Core_UI.getColor(1 : numel(sta_list), numel(sta_list)));
            ax(1) = subplot(3,1,1);
            
            if numel(sta_list) == 1
                Core_Utils.plotSep(p_time.getMatlabTime, pressure, '.-', 'Color', Core_UI.getColor(1));
            else
                Core_Utils.plotSep(p_time.getMatlabTime, pressure, '.-');
            end
            ylim( minMax(pressure) + [-2.5 2.5]);
            h = ylabel('Pres. [mbar]'); h.FontWeight = 'bold';

            outm = {};
            for r = 1 : numel(sta_list)
                outm{r} = ['P. ' sta_list(r).getMarkerName4Ch()];
            end
            [~, icons] = legend(outm, 'Location', 'NorthEastOutside', 'interpreter', 'none');
            n_entry = numel(outm);
            icons = icons(n_entry + 2 : 2 : end);
            for i = 1 : numel(icons)
                icons(i).MarkerSize = 18;
            end

            ax(2) = subplot(3,1,2);
            if numel(sta_list) == 1
                Core_Utils.plotSep(p_time.getMatlabTime, temperature, '.-', 'Color', Core_UI.getColor(2));
            else
                Core_Utils.plotSep(p_time.getMatlabTime, temperature, '.-');
            end
            ylim( minMax(temperature) + [-2.5 2.5]);
            h = ylabel('Temp. [deg C]'); h.FontWeight = 'bold';
            
            outm = {};
            for r = 1 : numel(sta_list)
                outm{r} = ['T. ' sta_list(r).getMarkerName4Ch()];
            end
            [~, icons] = legend(outm, 'Location', 'NorthEastOutside', 'interpreter', 'none');
            n_entry = numel(outm);
            icons = icons(n_entry + 2 : 2 : end);
            for i = 1 : numel(icons)
                icons(i).MarkerSize = 18;
            end
            
            ax(3) = subplot(3,1,3);
            if numel(sta_list) == 1
                Core_Utils.plotSep(p_time.getMatlabTime, humidity, '.-', 'Color', Core_UI.getColor(3));
            else
                Core_Utils.plotSep(p_time.getMatlabTime, humidity, '.-');
            end
            ylim([0 105]);
            h = ylabel('Hum. [%]'); h.FontWeight = 'bold';

            outm = {};
            for r = 1 : numel(sta_list)
                outm{r} = ['H. ' sta_list(r).getMarkerName4Ch()];
            end
            [~, icons] = legend(outm, 'Location', 'NorthEastOutside', 'interpreter', 'none');
            n_entry = numel(outm);
            icons = icons(n_entry + 2 : 2 : end);
            for i = 1 : numel(icons)
                icons(i).MarkerSize = 18;
            end

            linkaxes(ax, 'x');
            xlim([p_time.first.getMatlabTime() p_time.last.getMatlabTime()]);
            for i = 1:3
                grid(ax(i), 'minor');
                setTimeTicks(ax(i), 3, 'auto');                
            end
            
            Core_UI.beautifyFig(f);
            Core_UI.addExportMenu(f);
            Core_UI.addBeautifyMenu(f);
            f.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;            
        end

        function [tropo_az, tropo_modulus, gne, time_grad, time_grad_ref] = getTropoGradients(sta_list, max_n_val)
            % Get regularly sampled tropospheric gradients from out
            %
            % INPUT
            %   max_n_val      number of gradients to extract
            %
            % OUTPUT
            %   az             azimuth of the tropospheric gradient from North clockwise (deg)
            %   modulus        modulus of the tropospheric gradient [m]
            %   gne            North East gradients
            %   time_grad      GPS_Time of the exported gradients
            %   time_grad_ref  Nominal epoch of the searched gradient          
            %
            % SYNATAX
            %   [tropo_az, tropo_modulus, gne, time_grad, time_grad_ref] = GNSS_Station.getTropoGradients(sta_list, <max_n_val>);

            
            if nargin < 2
                max_n_val = 100;
            end
            
            % Get all the times
            time = GNSS_Station.getSyncTimeExpanded(sta_list, 30);
            time = time.getNominalTime;
            [t_ref, t0] = time.getRefTime;
            t_rate = median(diff(t_ref));
            arrow_rate = round(((t_ref(end) - t_ref(1)) / t_rate) / (max_n_val - 1));
            time_grad_ref = (0 : arrow_rate * t_rate : arrow_rate * t_rate * (max_n_val - 1))';

            % Supposing already synced gradients
            [tropo_ge, ~] = sta_list.getTropoPar('ge');
            [tropo_gn, t] = sta_list.getTropoPar('gn');
            
            if ~iscell(t)
                t = {t};
                tropo_gn = {tropo_gn};
                tropo_ge = {tropo_ge};
            end
            tropo_az = {};
            tropo_modulus = {};
            time_grad = {};
            gne = {};
            for r = 1 : numel(sta_list)
                time = t{r}.getNominalTime.getRefTime(double(t0));
                
                time_grad{r} = [];
                gne{r} = nan(numel(time_grad_ref), 2);
                for e = 1 : numel(time_grad_ref)
                    [~, id_min] = min(abs(time - time_grad_ref(e)));
                    if ~isempty(id_min)
                        time_grad{r} = [time_grad{r} time(id_min)];
                        gne{r}(e, 1) = tropo_gn{r}(id_min);
                        gne{r}(e, 2) = tropo_ge{r}(id_min);
                    end
                end
                
                time_grad{r} = time_grad{r}/86400 + double(t0);
                tropo_modulus{r} = double(hypot(gne{r}(:,1), gne{r}(:,2)));                
                tropo_az{r} = 90 - atan2d(gne{r}(:,1), gne{r}(:,2));
            end
            time_grad_ref = time_grad_ref/86400 + double(t0);
        end

        function fh_list = showTropoGradientsMR(sta_list)
            % Show arrow plot for gradients
            %
            % SYNATAX
            %   sta_list.showTropoGradientsMR()

            %% DEBUG sta_list = rec(~rec.isEmpty_mr);
            sta_list = sta_list(~sta_list.isEmpty_mr);
            
            if numel(sta_list) == 1
                fh_list = sta_list.showZtdAndGradients();
            else
                % Order Station by longitude
                enu = sta_list.getMedianPosENU;
                [~, id_sort] = sort(enu(:,2));
                sta_list = sta_list(id_sort);
                
                max_n_arrows = 100;
                
                [tropo_az, tropo_modulus, gne, time_grad, time_grad_ref] = sta_list.getTropoGradients(max_n_arrows);
                
                % Scale factor gradients
                g_scale = 1e3; % m -> mm
                
                f = figure('Visible', 'off');
                f.Name = sprintf('%03d: Gradients', f.Number); f.NumberTitle = 'off';
                fh_list = f;
                fig_name = 'All_Tropo_Gradients';
                f.UserData = struct('fig_name', fig_name);
                
                ax = axes;
                
                xlim([time_grad_ref(1) time_grad_ref(end)] + median(diff(time_grad_ref), 'omitnan') .* [-1 1]);
                ylim([0 numel(sta_list) + 1]);
                Core_UI.beautifyFig(f);
                drawnow
                
                plot_scale_factor = ax.PlotBoxAspectRatio(2) * diff(xlim) / diff(ylim);
                
                % Compute the 0.8 percentile modulus
                mod_perc = g_scale * perc(cell2mat(tropo_modulus'), 0.9);
                max_size = min(1, (2 * diff(xlim) / max_n_arrows) / plot_scale_factor);
                for r = 1 : numel(sta_list)
                    triPlot(time_grad{r}, ...
                        r * ones(numel(tropo_az{r}), 1), ...
                        tropo_az{r}, ...
                        tropo_modulus{r} * g_scale, ...
                        max_size, ...
                        max_size/mod_perc);
                end
                setTimeTicks();
                yticks(1 : numel(sta_list));
                yticklabels(sta_list.getMarkerName4Ch());
                cb = colorbar;
                title(sprintf('Tropospheric gradients [mm]\nNorth on top\\fontsize{5} \n'), 'FontName', 'Open Sans');
                Core_UI.beautifyFig(f);
                
                Core_UI.addExportMenu(f);
                Core_UI.addBeautifyMenu(f);
                f.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                
                colormap(flipud(Cmap.get('RdBu')));
            end
        end
        
        function fh_list = showTropoPar(sta_list, par_name, new_fig, sub_plot_nsat, flag_od)
            % Show a unique plot for all the stations given a certain data parameter
            %
            % SYNTAX
            %   fh_list = sta_list.showTropoPar(par_name, flag_new_fig, sub_plot_nsat, flag_od)
            
            % one function to rule them all

            fh_list = [];
            if nargin < 5 || isempty(flag_od) || flag_od
                [tropo, t, id_ko] = sta_list.getTropoPar(par_name);                
            else
                [tropo, t] = sta_list.getTropoPar(par_name);
                id_ko = [];
            end
            
            if ~iscell(tropo)
                tropo = {tropo};
                t = {t};
            end

            rec_ok = false(numel(sta_list), 1);
            for r = 1 : size(sta_list, 2)
                rec_ok(r) = ~isempty(tropo{r}) & any(~isnan(tropo{r}));
            end

            sta_list = sta_list(rec_ok);
            tropo = tropo(rec_ok);
            t = t(rec_ok);
            if ~isempty(id_ko)
                id_ko = id_ko(rec_ok);
            end
            
            flag_ok = false; for i = 1 : numel(tropo); flag_ok = flag_ok || any(tropo{i}); end
            
            if numel(sta_list) == 0 || ~flag_ok
                log = Core.getLogger();
                log.addError('The requested atmospheric parameter have not been computed or stored');
            else
                if nargin < 3
                    new_fig = true;
                end

                if nargin < 4 || isempty(sub_plot_nsat)
                    sub_plot_nsat = false;
                end
                
                nsat_is_empty = false;
                tmp = sta_list.getNumSat; 
                if iscell(tmp)
                    for i = 1 : numel(tmp)
                        nsat_is_empty = nsat_is_empty || ~any(tmp{i});
                    end
                else
                    nsat_is_empty = ~any(tmp);
                end
                sub_plot_nsat = sub_plot_nsat && ~nsat_is_empty;
                
                if isempty(tropo)
                    sta_list(1).out.log.addWarning([par_name ' and slants have not been computed']);
                else
                    tlim = [inf -inf];
                    if strcmpi(par_name, 'zwd')
                        dlim = [-1 50];
                    else
                        dlim = [inf -inf];
                    end
                    
                    if new_fig
                        cc = Core.getState.getConstellationCollector;
                        drawnow;
                        fh = figure('Visible', 'off'); fh.Name = sprintf('%03d: %s %s', fh.Number, par_name, cc.sys_c); fh.NumberTitle = 'off';
                        Core_UI.beautifyFig(fh);
                        drawnow;                        
                    else
                        fh = gcf;
                    end
                    
                    fh_list = [fh_list; fh];
                    if numel(sta_list) == 1 
                        % If I have only one receiver use as name the name of the receiver
                        fig_name = sprintf('%s_%s_%s', upper(par_name), sta_list.getMarkerName4Ch, sta_list.getTime.first.toString('yyyymmdd_HHMM'));
                    else
                        % If I have more than one receiver use as name the name of the project
                        fig_name = sprintf('%s_%s', upper(par_name), strrep(Core.getState.getPrjName,' ', '_'));
                    end
                    fh.UserData = struct('fig_name', fig_name);
                    
                    if new_fig
                        if sub_plot_nsat
                            ax1 = subplot(3,1,1:2);
                        else
                            ax1 = subplot(3,1,1:3);
                        end
                    else
                        try
                            ax1 = fh.Children(end);
                            subplot(ax1);
                        catch
                            ax1 = axes();
                        end
                    end
                    
                    if new_fig
                        old_legend = {};
                    else
                        l = legend;
                        old_legend = get(l,'String');
                        fh = gcf();
                    end
                    
                    drawnow;
                    e = 0;
                    set(0, 'CurrentFigure', fh);
                    for r = 1 : numel(sta_list)
                        rec = sta_list(r);
                        data_tmp = tropo{r};
                        if ~isempty(id_ko)
                            id_ko_tmp = id_ko{r};
                        else
                            id_ko_tmp = false(size(data_tmp));
                        end
                        mode = '.-';
                        if new_fig
                            if strcmpi(par_name, 'nsat')
                                Core_Utils.plotSep(t{r}.getMatlabTime(), zero2nan(data_tmp'), '.-', 'LineWidth', 2, 'Color', Core_UI.getColor(r, size(sta_list, 2))); hold on;
                            else
                                if any(id_ko_tmp)
                                    Core_Utils.plotSep(t{r}.getEpoch(find(id_ko_tmp)).getMatlabTime(), zero2nan(data_tmp(id_ko_tmp)').*1e2, mode, 'LineWidth', 2, 'Color', [0.9 0.9 0.9]); hold on;
                                    Core_Utils.plotSep(t{r}.getEpoch(find(~id_ko_tmp)).getMatlabTime(), zero2nan(data_tmp(~id_ko_tmp)').*1e2, mode, 'LineWidth', 2, 'Color', Core_UI.getColor(r, size(sta_list, 2))); hold on;
                                else
                                    Core_Utils.plotSep(t{r}.getMatlabTime(), zero2nan(data_tmp').*1e2, mode, 'LineWidth', 2, 'Color', Core_UI.getColor(r, size(sta_list, 2))); hold on;
                                end
                            end
                        else
                            if strcmpi(par_name, 'nsat')
                                plot(t{r}.getMatlabTime(), zero2nan(data_tmp'), '.-', 'LineWidth', 2); hold on;
                            else
                                if any(id_ko_tmp)
                                    Core_Utils.plotSep(t{r}.getEpoch(find(id_ko_tmp)).getMatlabTime(), zero2nan(data_tmp(id_ko_tmp)').*1e2, mode, 'LineWidth', 2, 'Color', [0.9 0.9 0.9]); hold on;
                                    Core_Utils.plotSep(t{r}.getEpoch(find(~id_ko_tmp)).getMatlabTime(), zero2nan(data_tmp(~id_ko_tmp)').*1e2, mode, 'LineWidth', 2); hold on;
                                else
                                    Core_Utils.plotSep(t{r}.getMatlabTime(), zero2nan(data_tmp').*1e2, mode, 'LineWidth', 2); hold on;
                                end
                            end
                        end
                        childs = findall(ax1.Children, 'type', 'line');
                        for c = 1 : length(childs)
                            tlim(1) = min([tlim(1), childs(c).XData]);
                            tlim(2) = max([tlim(2), childs(c).XData]);
                            
                            dlim(1) = min([dlim(1), childs(c).YData]);
                            dlim(2) = max([dlim(2), childs(c).YData]);
                        end
                        if strcmp(par_name, 'nsat') || ~any(id_ko_tmp)
                            e = e + 1;
                            outm{e} = rec(1).getMarkerName4Ch();
                        else
                            e = e + 1;
                            outm{e} = [rec(1).getMarkerName4Ch() ' outliers'];
                            e = e + 1;
                            outm{e} = rec(1).getMarkerName4Ch();
                        end
                    end
                    dspan = dlim(2) - dlim(1);
                    if dlim(1) < 0
                        dlim(1) = dlim(1) - 0.03 *dspan;
                    else
                        dlim(1) = max(0, dlim(1) - 0.03 *dspan);
                    end
                    dlim(2) = dlim(2) + 0.03 * dspan;
                    try
                        xlim(tlim);
                        ylim(dlim);
                    catch ex
                        % Grrr MATLAB is out of sync
                        Core_Utils.printEx(ex);
                        delete(fh);
                        keyboard;
                    end

                    outm = [old_legend, outm];
                    n_entry = numel(outm);

                    if n_entry > 50
                        loc = 'NorthEast';
                    else
                        loc = 'EastOutside';
                    end
                    [lh, icons] = legend(outm, 'Location', loc, 'interpreter', 'none');
                    icons = icons(n_entry + 2 : 2 : end);
                    
                    for i = 1 : numel(icons)
                        icons(i).MarkerSize = 18;
                    end
                    if n_entry > 50
                        lh.Visible = 'off';
                    end
                    

                    setTimeTicks(4);
                    h = ylabel([par_name ' [cm]']); h.FontWeight = 'bold';
                    grid on;
                    
                    % Generate Name
                    switch lower(par_name)
                        case 'ztd'
                            ttl = 'Estimated Zenith Total Delay (ZTD)';
                        case 'zwd'
                            ttl = 'Estimated Zenith Wet Delay (ZWD)';
                        case 'gn'
                            ttl = 'Estimated Wet Delay North Gradient';
                        case 'ge'
                            ttl = 'Estimated Wet Delay East Gradient';
                        case 'pwv'
                            ttl = 'Estimated Precipitable Water Vapour (PWV)';
                        case 'zhd'
                            ttl = 'Estimated Zenith Hydrostatic Delay (ZHD)';
                        case 'nsat'
                            ttl = 'Number of used satellites';
                    end
                    
                    h = title(ttl); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                    if sub_plot_nsat
                        ax2 = subplot(3,1,3);
                        for r = 1 : numel(sta_list)
                            rec = sta_list(r);
                            if new_fig
                                Core_Utils.plotSep(t{r}.getMatlabTime(), zero2nan(rec.getNumSat'), '.-', 'LineWidth', 2, 'Color', Core_UI.getColor(r, size(sta_list, 2))); hold on;
                            else
                                Core_Utils.plotSep(t{r}.getMatlabTime(), zero2nan(rec.getNumSat'), '.-', 'LineWidth', 2); hold on;
                            end
                        end
                        childs = ax2.Children;
                        tlim = [inf -inf];
                        dlim = [inf -inf];
                        for c = 1 : length(childs)
                            tlim(1) = min([tlim(1), childs(c).XData]);
                            tlim(2) = max([tlim(2), childs(c).XData]);
                            
                            dlim(1) = min([dlim(1), childs(c).YData]);
                            dlim(2) = max([dlim(2), childs(c).YData]);
                        end
                        dspan = dlim(2) - dlim(1);
                        if dlim(1) < 0
                            dlim(1) = dlim(1) - 1;
                        else
                            dlim(1) = max(0, dlim(1) - 1);
                        end
                        dlim(2) = dlim(2) + 1;
                        xlim(tlim);
                        ylim(dlim);
                        setTimeTicks(4);
                        if new_fig
                            h = ylabel(['# sat']); h.FontWeight = 'bold';
                            grid on;
                            linkaxes([ax1 ax2], 'x');
                        end                        
                        [lh, icons] = legend(outm, 'Location', loc, 'interpreter', 'none');
                        lh.Visible = 'off';
                    end
                end
                Core_UI.beautifyFig(fh);
                Core_UI.addExportMenu(fh);
                Core_UI.addBeautifyMenu(fh);
                fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
            end
        end

        function fh_list = showZwdSync(sta_list)
            % Show A matrix of all the ZWD
            % (If stations are more than 40 the names will be overposed)
            %
            % SYNTAX
            %   fh_list = sta_list.showZwdSync()
            
            [zwd, time] = sta_list.getZwdRes_mr;
            
            f = figure('Visible', 'off');
           
            fh_list = f;
            fig_name = sprintf('RecZwdSync');
            f.UserData = struct('fig_name', fig_name);
            
            imagesc(time.getMatlabTime, 1:size(zwd,2), zwd' * 1e2);
            colormap(Cmap.get('RdBu'));
            colorbar;
            yticks(1:size(zwd,2));
            yticklabels([sta_list.getMarkerName4Ch]);
            setTimeTicks();
            title(sprintf('ZWD [cm]\\fontsize{5} \n'));
            xlabel('Epochs');
            Core_UI.beautifyFig(f,'light');
                        
            Core_UI.addExportMenu(f);
            Core_UI.addBeautifyMenu(f);
            f.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
        end
        
        function fh_list = showZwdProcStatus(sta_list)
            % Show A matrix of all the ZWD
            % (If stations are more than 40 the names will be overposed)
            %
            % SYNTAX
            %   fh_list = sta_list.showZwdSync()
            
            [zwd, time] = sta_list.getZwdRes_mr;
            
            f = figure('Visible', 'off');
            av_res = ~isnan(zwd);
            av_data = false(size(av_res));
           
            fr = Core.getRinLists();
            time_matlab = time.getMatlabTime;
            rate = time.getRate /86400; %rate in days
            for r = 1 : length(sta_list)
                for fl = 1 : fr(r).first_epoch.length
                    idx_epoch = time_matlab > (fr(r).first_epoch.getEpoch(fl).getMatlabTime -rate +eps) & time_matlab < (fr(r).last_epoch.getEpoch(fl).getMatlabTime +rate -eps);
                    av_data(idx_epoch,r) = true;
                end
            end
                
            plot_data = zeros(size(av_data));
            plot_data(av_data) = 1;
            plot_data(av_res) = 2;
            fh_list = f;
            fig_name = sprintf('Availability stat');
            f.UserData = struct('fig_name', fig_name);
            
            imagesc(time.getMatlabTime, 1:size(zwd,2), plot_data');
            %colormap(Cmap.get('RdBu'));
            colormap(Cmap.getColor(1:3,3, 'Pastel1'))
            caxis([0 3]); cb = colorbar;
            cb.Ticks = [0.5 1.5 2.5]; cb.TickLabels = {'No data', 'No solution', 'ZWD Ok'};
            if verLessThan( 'MATLAB', '9.5' )
                set(gca,'YTick',1:size(zwd,2));
                set(gca,'YTickLabel',[sta_list.getMarkerName4Ch]);
            else
                yticks(1:size(zwd,2));
                yticklabels([sta_list.getMarkerName4Ch]);
            end
            setTimeTicks();
            title(sprintf('ZWD processing status\\fontsize{5} \n'));
            xlabel('Epochs');
            Core_UI.beautifyFig(f,'light');
                        
            Core_UI.addExportMenu(f);
            Core_UI.addBeautifyMenu(f);
            f.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
        end
        
        function fh_list = showNSat(sta_list, new_fig)
            % Show total number of satellites in view (epoch by epoch) for each satellite
            %
            % SYNTAX:
            %   fh_list = sta_list.showNSat()
            if nargin == 1
                new_fig = true;
            end
            fh_list = sta_list.showTropoPar('nsat', new_fig, false);
        end

        function fh_list = showNSatSS(sta_list, flag_smooth)
            % Show total number of satellites in view (epoch by epoch) for each satellite
            %
            % SYNTAX:
            %   fh_list = sta_list.showNSatSS()

            fh_list = [];
            sta_list = sta_list(~(sta_list.isEmpty_mr));            
            for r = 1 : numel(sta_list)
                if nargin == 2
                    if ~(isempty(sta_list(r).out) || sta_list(r).out.isEmpty)
                        fh_list = [fh_list; sta_list(r).out.showNSatSS(flag_smooth)];
                    else
                        fh_list = [fh_list; sta_list(r).work.showNSatSS(flag_smooth)];
                    end
                else
                    if ~(isempty(sta_list(r).out) || sta_list(r).out.isEmpty)
                        fh_list = [fh_list; sta_list(r).out.showNSatSS()];
                    else
                        fh_list = [fh_list; sta_list(r).work.showNSatSS()];
                    end
                end
            end
        end

        function fh_list = showZhd(sta_list, new_fig, sub_plot_nsat, flag_od)
            % Display ZHD values
            %
            % INPUT:
            %   new_fig         flag to specify to open a new figure (default = true)
            %   sub_plot_nsat   flag to specify to subplot #sat      (default = true)
            %   flag_od         flag to disable outlier detection
            %
            % SYNTAX:
            %   sta_list.showZhd(<new_fig = true>, <sub_plot_nsat = true>, <flag_od = false>)

            if nargin <= 1 || isempty(new_fig)
                new_fig = true;
            end
            if nargin <= 2 || isempty(sub_plot_nsat)
                sub_plot_nsat = true;
            end
            if nargin <= 3 || isempty(flag_od)
                flag_od = false;
            end
            
            fh_list = sta_list.showTropoPar('ZHD', new_fig, sub_plot_nsat, flag_od);
        end

        function fh_list = showZwd(sta_list, new_fig, sub_plot_nsat, flag_od)
            % Display ZWD values
            %
            % INPUT:
            %   new_fig         flag to specify to open a new figure (default = true)
            %   sub_plot_nsat   flag to specify to subplot #sat      (default = true)
            %   flag_od         flag to disable outlier detection
            %
            % SYNTAX:
            %   sta_list.showZwd(<new_fig = true>, <sub_plot_nsat = true>, <flag_od = false>)

            if nargin <= 1 || isempty(new_fig)
                new_fig = true;
            end
            if nargin <= 2 || isempty(sub_plot_nsat)
                sub_plot_nsat = true;
            end
            if nargin <= 3 || isempty(flag_od)
                flag_od = false;
            end
            
            fh_list = sta_list.showTropoPar('ZWD', new_fig, sub_plot_nsat, flag_od);
        end

        function fh_list = showPwv(sta_list, new_fig, sub_plot_nsat, flag_od)
            % Display PWV values
            %
            % INPUT:
            %   new_fig         flag to specify to open a new figure (default = true)
            %   sub_plot_nsat   flag to specify to subplot #sat      (default = true)
            %   flag_od         flag to disable outlier detection
            %
            % SYNTAX:
            %   sta_list.showPwv(<new_fig = true>, <sub_plot_nsat = true>, <flag_od = false>)

            if nargin <= 1 || isempty(new_fig)
                new_fig = true;
            end
            if nargin <= 2 || isempty(sub_plot_nsat)
                sub_plot_nsat = true;
            end
            if nargin <= 3 || isempty(flag_od)
                flag_od = false;
            end
            
            fh_list = sta_list.showTropoPar('PWV', new_fig, sub_plot_nsat, flag_od);
        end

        function fh_list = showZtd(sta_list, new_fig, sub_plot_nsat, flag_od)
            % Display ZTD values
            %
            % INPUT:
            %   new_fig         flag to specify to open a new figure (default = true)
            %   sub_plot_nsat   flag to specify to subplot #sat      (default = true)
            %   flag_od         flag to disable outlier detection
            %
            % SYNTAX:
            %   sta_list.showZtd(<new_fig = true>, <sub_plot_nsat = true>, <flag_od = false>)

            if nargin <= 1 || isempty(new_fig)
                new_fig = true;
            end
            if nargin <= 2 || isempty(sub_plot_nsat)
                sub_plot_nsat = true;
            end
            if nargin <= 3 || isempty(flag_od)
                flag_od = false;
            end
            
            fh_list = sta_list.showTropoPar('ZTD', new_fig, sub_plot_nsat, flag_od);
        end

        function fh_list = showZtdAndGradients(sta_list, new_fig, sub_plot_nsat, flag_od)
            % Display ZTD values
            % For each receiver use a single plot and overimpose the gradients
            %
            %
            % INPUT:
            %   new_fig         flag to specify to open a new figure (default = true)
            %   sub_plot_nsat   flag to specify to subplot #sat      (default = true)
            %   flag_od         flag to disable outlier detection
            %
            % SYNTAX:
            %   sta_list.showZtdAndGradients(<new_fig = true>, <sub_plot_nsat = true>, <flag_od = false>)

            if nargin <= 1 || isempty(new_fig)
                new_fig = true;
            end
            if nargin <= 2 || isempty(sub_plot_nsat)
                sub_plot_nsat = true;
            end
            if nargin <= 3 || isempty(flag_od)
                flag_od = false;
            end
            
            sta_list = sta_list(~sta_list.isEmpty_mr);

            g_scale = 1e2;
            max_n_arrows = 120;
            fh_list = sta_list.showTropoPar('ZTD', new_fig, sub_plot_nsat, flag_od);
            fh_list.UserData.fig_name = [fh_list.UserData.fig_name 'AndGrad'];
            for r = 1 : numel(sta_list)                                
                [ztd, time] = sta_list(r).getTropoPar('ZTD');
                ax = fh_list(end).Children(end);
                subplot(ax);
                
                [tropo_az, tropo_modulus, gne, time_grad, ~] = sta_list(r).getTropoGradients(max_n_arrows);
                ztd = interp1(time.getMatlabTime, ztd, time_grad{1});
                
                plot_scale_factor = ax.PlotBoxAspectRatio(2) * diff(xlim) / diff(ylim);
                mod_perc = g_scale * perc(cell2mat(tropo_modulus'), 0.9);
                max_size = (3 * diff(xlim) / max_n_arrows) / plot_scale_factor;
                triPlot(time_grad{1}, ...
                    ztd * g_scale, ...
                    tropo_az{1}, ...
                    tropo_modulus{1} * g_scale, ...
                    max_size);                
            end
        end
        
        function fh_list = showGn(sta_list, new_fig, sub_plot_nsat, flag_od)
            % Display ZTD Gradiet North values
            %
            % INPUT:
            %   new_fig         flag to specify to open a new figure (default = true)
            %   sub_plot_nsat   flag to specify to subplot #sat      (default = true)
            %   flag_od         flag to disable outlier detection
            %
            % SYNTAX:
            %   sta_list.showGn(<new_fig = true>, <sub_plot_nsat = true>)

            if nargin <= 1 || isempty(new_fig)
                new_fig = true;
            end
            if nargin <= 2 || isempty(sub_plot_nsat)
                sub_plot_nsat = true;
            end
            if nargin <= 3 || isempty(flag_od)
                flag_od = false;
            end
            fh_list = sta_list.showTropoPar('GN', new_fig, sub_plot_nsat, flag_od);
        end

        function fh_list = showGe(sta_list, new_fig, sub_plot_nsat, flag_od)
            % Display ZTD Gradiet East values
            %
            % INPUT:
            %   new_fig         flag to specify to open a new figure (default = true)
            %   sub_plot_nsat   flag to specify to subplot #sat      (default = true)
            %   flag_od         flag to disable outlier detection
            %
            % SYNTAX:
            %   sta_list.showGe(<new_fig = true>, <sub_plot_nsat = true>)

            if nargin <= 1 || isempty(new_fig)
                new_fig = true;
            end
            if nargin <= 2 || isempty(sub_plot_nsat)
                sub_plot_nsat = true;
            end
            if nargin <= 2 || isempty(flag_od)
                flag_od = false;
            end
            fh_list = sta_list.showTropoPar('GE', new_fig, sub_plot_nsat, flag_od);
        end
        
        function fh_list = showZtdVsHeight(sta_list, degree)
            % Show Median ZTD of n_receivers vs Orthometric height
            %
            % SYNTAX
            %   sta_list.showZtdVsHeight();
            f = figure('Visible', 'off');
            
            fh_list = f;
            fig_name = sprintf('ZtdVsHeight');
            f.UserData = struct('fig_name', fig_name);
                        
            med_ztd = median(sta_list.getZtd_mr * 1e2, 'omitnan')';
            ax1 = subplot(2,1,1);
            [lat, lon, ~, h_o] = Coordinates.fromXYZ(sta_list.getMedianPosXYZ()).getGeodetic;
            [~, id_sort] = sort(lat);
            plot(h_o(id_sort), med_ztd(id_sort), '.', 'MarkerSize', 20); hold on;
            %scatter(h_o(id_sort), med_ztd(id_sort), 200, lat(id_sort), '.'); 
            %id=(lat(id_sort)/pi*180)>41; plot(h_o(id_sort(id)), med_ztd(id_sort(id)), '.', 'MarkerSize', 20); hold on;
            %plot(h_o(id_sort(id)), med_ztd(id_sort(id)), '.', 'MarkerSize', 20); hold on;
            colormap(Cmap.get('viridis', numel(lat)));
            hold on;
            ylabel('Median ZTD [cm]');
            xlabel('Elevation [m]');
            title('ZTD vs Height')

            if nargin == 1
                degree = 2;
            end
            y_out = Core_Utils.interp1LS(h_o, med_ztd, degree, h_o);
            plot(sort(h_o), Core_Utils.interp1LS(h_o, med_ztd, degree, sort(h_o)), '-', 'Color', Core_UI.COLOR_ORDER(3,:), 'LineWidth', 2);
            grid on
            ax1.FontSize = 16;
            
            ax2 = subplot(2,1,2);
            plot(h_o, med_ztd - y_out, '.', 'MarkerSize', 20); hold on
            %plot(h_o(id_sort(id)), med_ztd(id_sort(id)) - y_out(id_sort(id)), '.', 'MarkerSize', 20); hold on;

            ylabel('Residual [cm]');
            xlabel('Elevation [m]');
            title('Reduced ZTD vs Height', 'FontName', 'Open Sans')

            sta_strange = find(abs(med_ztd - y_out) > 8);
            if ~isempty(sta_strange)
                Core.getLogger.addMessage('Strange station detected');
                for s = 1 : numel(sta_strange)
                    Core.getLogger.addMessage(sprintf(' %d - %s', sta_strange(s), sta_list(sta_strange(s)).getMarkerName()));
                end
            end
            grid on
            ax2.FontSize = 16;
            Core_UI.beautifyFig(gcf);
            Core_UI.addExportMenu(gcf);             
            Core_UI.addBeautifyMenu(gcf);             
            f.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
        end

        function fh_list = showZwdVsHeight(sta_list, degree)
            % Show Median ZTD of n_receivers vs Orthometric height
            %
            % SYNTAX
            %   sta_list.showZwdVsHeight();
            f = figure('Visible', 'off');
            
            fh_list = f;
            fig_name = sprintf('ZwdVsHeight');
            f.UserData = struct('fig_name', fig_name);
            
            med_zwd = median(sta_list.getZwd_mr * 1e2, 'omitnan')';
            
            ax1 = subplot(2,1,1);
            [~, ~, ~, h_o] = Coordinates.fromXYZ(sta_list.getMedianPosXYZ()).getGeodetic;
            plot(h_o, med_zwd, '.', 'MarkerSize', 20); hold on;
            ylabel('Median ZWD [cm]');
            xlabel('Elevation [m]');
            title('ZWD vs Height')

            if nargin == 1
                degree = 2;
            end
            y_out = Core_Utils.interp1LS(h_o, med_zwd, degree, h_o);
            plot(sort(h_o), Core_Utils.interp1LS(h_o, med_zwd, degree, sort(h_o)), '-', 'Color', Core_UI.COLOR_ORDER(3,:), 'LineWidth', 2);
            grid on
            ax1.FontSize = 16;
            
            ax2 = subplot(2,1,2);
            plot(h_o, med_zwd - y_out, '.', 'MarkerSize', 20);

            ylabel('Residual [cm]');
            xlabel('Elevation [m]');
            title('Reduced ZWD vs Height', 'FontName', 'Open Sans')

            sta_strange = find(abs(med_zwd - y_out) > 8);
            if ~isempty(sta_strange)
                Core.getLogger.addMessage('Strange station detected');
                for s = 1 : numel(sta_strange)
                    Core.getLogger.addMessage(sprintf(' %d - %s', sta_strange(s), sta_list(sta_strange(s)).getMarkerName()));
                end
            end
            grid on
            Core_UI.beautifyFig(gcf);
            Core_UI.addExportMenu(gcf);             
            Core_UI.addBeautifyMenu(gcf);             
            f.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
        end

        function fh_list = showMedianTropoPar(this, par_name, new_fig)
            % one function to rule them all
            
            fh_list = [];
            
            rec_ok = false(size(this,2), 1);
            for r = 1 : size(this, 2)
                rec_ok(r) = any(~isnan(this(:,r).out.getZtd));
            end
            sta_list = this(:, rec_ok);

            if nargin < 3
                new_fig = true;
            end

            switch lower(par_name)
                case 'ztd'
                    [tropo] = sta_list(1).out.getZtd();
                case 'zwd'
                    [tropo] = sta_list(1).out.getZwd();
                case 'pwv'
                    [tropo] = sta_list(1).out.getPwv();
                case 'zhd'
                    [tropo] = sta_list(1).out.getAprZhd();
            end

            if ~iscell(tropo)
                tropo = {tropo};
            end
            if isempty(tropo)
                sta_list(1).out.log.addWarning([par_name ' and slants have not been computed']);
            else
                if new_fig
                    f = figure('Visible', 'off'); f.Name = sprintf('%03d: Median %s %s', f.Number, par_name, sta_list(1).out.getCC.sys_c); f.NumberTitle = 'off';
                    old_legend = {};
                else
                    f = gcf;
                    l = legend;
                    old_legend = get(l,'String');
                end
                
                fh_list = f;
                % If I have more than one receiver use as name the name of the project
                fig_name = sprintf('%s_Median_%s', upper(par_name), strrep(Core.getState.getPrjName,' ', '_'));
                f.UserData = struct('fig_name', fig_name);
                
                for r = 1 : numel(sta_list)
                    rec = sta_list(~sta_list(r).isEmpty, r);
                    if ~isempty(rec)
                        switch lower(par_name)
                            case 'ztd'
                                [tropo] = rec.out.getZtd();
                            case 'zwd'
                                [tropo] = rec.out.getZwd();
                            case 'pwv'
                                [tropo] = rec.out.getPwv();
                            case 'zhd'
                                [tropo] = rec.out.getAprZhd();
                        end
                        [~, ~, ~, h_o] = rec(1).out.getPosGeodetic();
                        if new_fig
                            plot(h_o, zero2nan(median(tropo,'omitnan')), '.', 'MarkerSize', 35, 'LineWidth', 4, 'Color', Core_UI.getColor(r, size(sta_list, 2))); hold on;
                        else
                            plot(h_o, zero2nan(median(tropo,'omitnan')), '.', 'MarkerSize', 35, 'LineWidth', 4); hold on;
                        end
                        outm{r} = rec(1).getMarkerName();
                        h_ortho(r) = h_o;
                        med_tropo(r) = median(tropo,'omitnan');
                    else
                        h_ortho(r) = nan;
                        med_tropo(r) = nan;
                    end
                end

                h_ortho(med_tropo == 0) = nan;
                med_tropo = zero2nan(med_tropo);
                degree = 2;
                h_grid = min(noNaN(h_ortho)) :  min(10, diff(minMax(noNaN(h_ortho)))/100) : max(noNaN(h_ortho));
                h_component = Core_Utils.interp1LS(noNaN(h_ortho), noNaN(med_tropo), degree, h_grid);
                plot(h_grid, h_component, '-k', 'LineWidth', 2);

                outm = [old_legend, outm];
                [~, icons] = legend(outm, 'Location', 'NorthEastOutside', 'interpreter', 'none');
                n_entry = numel(outm);
                icons = icons(n_entry + 2 : 2 : end);

                for i = 1 : numel(icons)
                    icons(i).MarkerSize = 18;
                end

                %ylim(yl);
                %xlim(t(time_start) + [0 win_size-1] ./ 86400);
                h = ylabel([par_name ' [m]']); h.FontWeight = 'bold';
                h = xlabel('Elevation [m]'); h.FontWeight = 'bold';
                grid on;
                h = title(['Median Receiver ' par_name]); h.FontWeight = 'bold'; %h.Units = 'pixels'; h.Position(2) = h.Position(2) + 8; h.Units = 'data';
                
                Core_UI.beautifyFig(gcf);
                Core_UI.addExportMenu(gcf);             
                Core_UI.addBeautifyMenu(gcf);             
                f.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
            end
        end

        function fh_list = showMedianZhd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showMedianTropoPar('ZHD', new_fig);
        end

        function fh_list = showMedianZwd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showMedianTropoPar('ZWD', new_fig);
        end

        function fh_list = showMedianZtd(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showMedianTropoPar('ZTD', new_fig);
        end

        function fh_list = showMedianPwv(this, new_fig)
            if nargin == 1
                new_fig = true;
            end
            fh_list = this.showMedianTropoPar('PWV', new_fig);
        end

        function fh_list = showZtdSlantRes_p(this, time_start, time_stop)
            fh_list = [];
            for r = 1 : size(this, 2)
                if nargin == 1
                    time_start = this(r).out.getTime.first;
                    time_stop = this(r).out.getTime.last;
                end
                ztd = this(r).out.getZtd();
                slant_td = this(r).out.getSlantTD();
                if isempty(ztd) || ~any(slant_td(:))
                    this.log.addWarning('ZTD and slants have not been computed');
                else
                    fh_list = [fh_list; this(r).out.showZtdSlantRes_p(time_start, time_stop)];
                end
            end
        end
        
        function fh_list = showBaselineENU(sta_list, baseline_ids, flag_add_coo, n_obs)
            % Function to plot baseline between 2 or more stations
            %
            % INPUT:
            %   sta_list                 list of GNSS_Station objects
            %   baseline_ids/ref_id      n_baseline x 2 - couple of id in sta_list to be used
            %   flag_add_coo             use external / internal / high rate cordinates
            %   n_obs                    use only the last n_obs
            %
            % SYNTAX
            %   sta_list.showBaselineENU(<baseline_ids = []>, <flag_add_coo>, <n_obs>)
            if nargin < 4 || isempty(n_obs) || isnan(n_obs)
                n_obs = 0;
            end
            if nargin < 3
                flag_add_coo = 0;
            end
            if flag_add_coo < 0
                tmp_id = 1:numel(sta_list); % find valid
            else
                tmp_id = find(~sta_list.isEmptyOut_mr); % find valid
            end
            
            if flag_add_coo >= 0
                sta_list = sta_list(~sta_list.isEmptyOut_mr);
            end
            out_list = [sta_list.out];
            
            id2valid = nan(numel(sta_list), 1);
            id2valid(tmp_id) = 1:numel(tmp_id);
            if nargin < 2
                baseline_ids = [ones(numel(out_list)-1,1) (2:numel(out_list))'];
            else
                baseline_ids(:) = id2valid(baseline_ids(:));
            end
            
            fh_list = out_list.showBaselineENU(baseline_ids, flag_add_coo, n_obs);
        end
        
        function fh_list = showBaselinePlanarUp(sta_list, baseline_ids, flag_add_coo, n_obs)
            % Function to plot baseline between 2 or more stations
            %
            % INPUT:
            %   sta_list                 list of GNSS_Station objects
            %   baseline_ids/ref_id      n_baseline x 2 - couple of id in sta_list to be used
            %   flag_add_coo             use external / internal / high rate cordinates
            %   n_obs                    use only the last n_obs
            %
            % SYNTAX
            %   sta_list.showBaselinePlanarUp(<baseline_ids = []>, <flag_add_coo>, <n_obs>)
            if flag_add_coo < 0
                tmp_id = 1:numel(sta_list); % find valid
            else
                tmp_id = find(~sta_list.isEmptyOut_mr); % find valid
            end
            id2valid = nan(numel(sta_list), 1);
            id2valid(tmp_id) = 1:numel(tmp_id);
            if flag_add_coo >= 0
                sta_list = sta_list(~sta_list.isEmptyOut_mr);
            end
            out_list = [sta_list.out];
            if nargin < 4 || isempty(n_obs) || isnan(n_obs)
                n_obs = 0;
            end
            if nargin < 3
                flag_add_coo = 0;
            end
            if nargin < 2
                baseline_ids = [ones(numel(out_list)-1,1) (2:numel(out_list))'];
            else
                baseline_ids(:) = id2valid(baseline_ids(:));
            end
            fh_list = out_list.showBaselinePlanarUp(baseline_ids, flag_add_coo, n_obs);
        end
                                
        function fh_list = showMapDtmWithCloseRaob(sta_list)
            [id_rds, lat, lon, appr_dist] = sta_list.getCloseRaobIdList(Radiosonde.MAX_DIST); % Get RAOB within 150Km
            fh_list = sta_list.showMapDtm([], [], lat, lon);
            fh.UserData.fig_name = 'RecRaobMapDtm';
            
            [x, y] = m_ll2xy(lon, lat);
            
            % Label BG (in background w.r.t. the point)
            for r = 1 : numel(id_rds)
                name = Radiosonde.getRaobName(id_rds(r));
                text(x(r), y(r), ['     ' name ' '], ...
                    'FontWeight', 'bold', 'FontSize', 12, 'Color', [1 1 1], ...
                    'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                    'Margin', 2, 'LineWidth', 2, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'middle');
            end
            
            rds = plot(x, y, 's', ...
                'MarkerSize', 19, ...
                'MarkerEdgeColor', [0 0 0], ...
                'MarkerFaceColor', [0 0 0]);
            rds = plot(x, y, 's', ...
                'MarkerSize', 17, ...
                'MarkerEdgeColor', [0 0 0], ...
                'MarkerFaceColor', Core_UI.ORANGE, ...
                'UserData', 'RAOB_point');
            plot(x(:), y(:), '.k', 'MarkerSize', 5);
            
            % Label
            for r = 1 : numel(id_rds)
                name = Radiosonde.getRaobName(id_rds(r));
                text(x(r), y(r), ['     ' name ' '], ...
                    'FontWeight', 'bold', 'FontSize', 12, 'Color', [0 0 0], ...                    
                    'Margin', 2, 'LineWidth', 2, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'middle');
            end
            
            Core_UI.beautifyFig(fh_list, 'dark');            
        end
        
        function fh_list = showMapGoogleWithCloseRaob(sta_list)
            [id_rds, lat, lon, appr_dist] = sta_list.getCloseRaobIdList(Radiosonde.MAX_DIST); % Get RAOB within 150Km
            fh_list = sta_list.showMapGoogle([], lat, lon);
            fh.UserData.fig_name = 'RecRaobMapDtm';            
            [x, y] = m_ll2xy(lon, lat);
            
            % Label BG (in background w.r.t. the point)
            for r = 1 : numel(id_rds)
                name = Radiosonde.getRaobName(id_rds(r));
                text(x(r), y(r), ['     ' name ' '], ...
                    'FontWeight', 'bold', 'FontSize', 12, 'Color', [1 1 1], ...
                    'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                    'Margin', 2, 'LineWidth', 2, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'middle');
            end
            
            rds = plot(x, y, 's', ...
                'MarkerSize', 19, ...
                'MarkerEdgeColor', [0 0 0], ...
                'MarkerFaceColor', [0 0 0]);
            rds = plot(x, y, 's', ...
                'MarkerSize', 17, ...
                'MarkerEdgeColor', [0 0 0], ...
                'MarkerFaceColor', Core_UI.ORANGE, ...
                'UserData', 'RAOB_point');
            plot(x(:), y(:), '.k', 'MarkerSize', 5);
            
            % Label
            for r = 1 : numel(id_rds)
                name = Radiosonde.getRaobName(id_rds(r));
                text(x(r), y(r), ['     ' name ' '], ...
                    'FontWeight', 'bold', 'FontSize', 12, 'Color', [0 0 0], ...
                    'Margin', 2, 'LineWidth', 2, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'middle');
            end
            
            Core_UI.beautifyFig(fh_list, 'dark');
        end        
        
        function [fh_list, m_diff, s_diff, rds] = showRaobZtdValidation(sta_list, rds_list)
            % Compute and show comparison with radiosondes from weather.uwyo.edu
            % given region list, and station id (as cell arrays)
            %
            % INPUT
            %   sta_list        list of gnss receivers
            %   rds_list        cell array of string containing the radiosonde ID as used at "http://weather.uwyo.edu/upperair/sounding.html"
            %                   if it is empty find the closest in the area
            %
            % OUTPUT
            %   fh_list         figures handles
            %   m_diff          mean of the ZTD differences
            %   s_diff          std of the ZTD differences
            %
            % SYNTAX
            %  [m_diff, s_diff] = sta_list.showRaobZtdValidation(station_list);
            %
            % EXAMPLE
            %  % testing geonet full network
            %  [m_diff, s_diff] = sta_list.showRaobZtdValidation();
            %  [m_diff, s_diff] = sta_list.showRaobZtdValidation(Radiosonde.JAPAN_STATION);
            %
            % SEE ALSO
            %   Radiosonde GNSS_Station.getRaobZtdValidation
            if nargin < 2
                rds_list = [];
            end
            [m_diff, s_diff, rds, fh_list] = sta_list.getRaobZtdValidation(rds_list, true);
        end
        
        function [m_diff, s_diff, rds, fh_list] = getRaobZtdValidation(sta_list, rds_list, flag_show)
            % Compute a comparison with radiosondes from weather.uwyo.edu
            % given region list, and station id (as cell arrays)
            %
            % INPUT
            %   sta_list        list of gnss receivers
            %   rds_list        cell array of string containing the radiosonde ID as used at "http://weather.uwyo.edu/upperair/sounding.html"
            %
            % OUTPUT
            %   m_diff          mean of the ZTD differences
            %   s_diff          std of the ZTD differences
            %
            % SYNTAX
            %  [m_diff, s_diff] = sta_list.getRaobZtdValidation(station_list);
            %
            % EXAMPLE
            %  % testing geonet full network
            %  [m_diff, s_diff] = sta_list.getRaobZtdValidation(Radiosonde.JAPAN_STATION);
            %
            % SEE ALSO
            %   Radiosonde
            
            log = Core.getLogger();
            fh_list = [];
            if nargin < 3
                flag_show = false;
            end
            
            if nargin < 2 || isempty(rds_list)
                rds_list = sta_list.getCloseRaobIdList(Radiosonde.MAX_DIST);
            end
            
            if isempty(rds_list)
                log.addError('No radiosonde found, Validation Failed');
            else
                log.addMarkedMessage('Retrieving data, please wait...');
                % Get time limits
                p_time = sta_list.getSyncTimeExpanded(sta_list);
                start_time = GPS_Time(floor(p_time.first.getMatlabTime * 2)/2);
                stop_time = GPS_Time(ceil(p_time.last.getMatlabTime * 2)/2);
                
                % Download radiosondes
                rds = Radiosonde.fromList(rds_list, start_time, stop_time);
                
                % Get closer GNSS station
                [id_rec, d3d, dup] = sta_list.getCloserRec(rds.getLat(), rds.getLon(), rds.getHeight());
                gnss_list = sta_list(id_rec);
                
                % Interpolate ZTD
                log.addMarkedMessage('Get GNSS interpolated ZTD @ radiosonde locations');
                if numel(sta_list) > 1
                    warning off;
                    [ztd, ztd_height_correction, time] = sta_list.getTropoInterp('ZTD', rds.getLat(), rds.getLon(), rds.getHeight(), 300);
                    warning on;
                else
                    log.addWarning('Interpolation is not possible with just one station!!!');
                    ztd_height_correction = 0;
                    [ztd, time] = sta_list(id_rec(1)).getZtd_mr();
                    ztd = ztd * 1e2;
                end
                
                %             [lat, ~, ~, h_o] = Coordinates.fromXYZ(sta_list.getMedianPosXYZ()).getGeodetic;
                %             [~, id_sort] = sort(lat);
                %             id_north = (lat(id_sort) / pi * 180) > 41;
                %             [ztd_n, ztd_height_correction_n] = sta_list(id_sort(id_north)).getTropoInterp('ZTD', rds.getLat(), rds.getLon(), rds.getHeight());
                %             [ztd_s, ztd_height_correction_s] = sta_list(id_sort(~id_north)).getTropoInterp('ZTD', rds.getLat(), rds.getLon(), rds.getHeight());
                %             ztd_ns(:,rds.getLat > 41) = ztd_n(:,rds.getLat > 41);
                %             ztd_ns(:,rds.getLat <= 41) = ztd_s(:,rds.getLat <= 41);
                
                % Compute values
                log.addMonoMessage(sprintf('\n---------------------------------------------------------------------------------------'));
                [m_diff, s_diff, m_diff_sta, s_diff_sta, m_diff_sta_hcorr, s_diff_sta_hcorr] = deal(nan(numel(rds), 1));
                log.addMonoMessage(sprintf(' ZTD Radiosonde Validation (vs interpolated ZTD)\n---------------------------------------------------------------------------------------'));
                log.addMonoMessage(sprintf('                                Closer              Elevation     '));
                log.addMonoMessage(sprintf('       Mean          Std         GNSS    Dist [km]  diff. [m]  Radiosonde Station'));
                log.addMonoMessage(sprintf('---------------------------------------------------------------------------------------'));
                for s = 1 : numel(rds)
                    %                 if rds(s).getLat > 41
                    %                     ztd = ztd_n;
                    %                     ztd_height_correction = ztd_height_correction_n;
                    %                 else
                    %                     ztd = ztd_s;
                    %                     ztd_height_correction = ztd_height_correction_s;
                    %                 end
                    
                    % radiosondes
                    [ztd_rds, time_rds] = rds(s).getZtd();
                    
                    id_min = zeros(time_rds.length, 1);
                    ztd_diff = nan(time_rds.length, 1);
                    for e = 1 : time_rds.length
                        [t_min, id_min(e)] = min(abs(time - time_rds.getEpoch(e)));
                        if t_min < 3600
                            ztd_diff(e) = ztd_rds(e) - (ztd(id_min(e),s) + ztd_height_correction(s));
                        end
                    end
                    
                    m_diff(s) = mean(ztd_diff, 1, 'omitnan');
                    s_diff(s) = std(ztd_diff, 1, 'omitnan');
                    log.addMonoMessage(sprintf('%2d)  %6.2f cm    %6.2f cm      %4s  %9.1f   %9.1f   "%s"', s, m_diff(s), s_diff(s), sta_list(id_rec(s)).getMarkerName4Ch, round(d3d(s) / 1e3), dup(s), rds(s).getName()));
                end
                log.addMonoMessage(sprintf('\n---------------------------------------------------------------------------------------'));
                if flag_show
                    
                    log.addMonoMessage(sprintf(' ZTD Radiosonde Validation (vs height corrected closest station)\n---------------------------------------------------------------------------------------'));
                    log.addMonoMessage(sprintf('                                Closer              Elevation     '));
                    log.addMonoMessage(sprintf('       Mean          Std         GNSS    Dist [km]  diff. [m]  Radiosonde Station\n'));
                    log.addMonoMessage(sprintf('---------------------------------------------------------------------------------------'));
                    % Plot comparisons
                    for s = 1 : numel(rds)
                        [s_ztd, s_time] = sta_list(id_rec(s)).getZtd_mr();
                        % radiosondes
                        [ztd_rds, time_rds] = rds(s).getZtd();
                        
                        fh = figure('Visible', 'off');
                        Core_UI.beautifyFig(fh);
                        fh.Name = sprintf('%03d: Rds %s', fh.Number, rds(s).sta_num); fh.NumberTitle = 'off';
                        fh_list = [fh_list; fh]; %#ok<AGROW>
                        fig_name = sprintf('Raob_ZTD_%s_validation', rds(s).sta_num);
                        fh.UserData = struct('fig_name', fig_name);
                        
                        % interpolated ZTD
                        Core_Utils.plotSep(time.getMatlabTime, ztd(:,s) + ztd_height_correction(s), '.-', 'LineWidth', 2);
                        hold on;
                        
                        % closer ZTD
                        Core_Utils.plotSep(s_time.getMatlabTime, s_ztd * 1e2, '.-', 'LineWidth', 2);
                        
                        % get Radiosonde coordinates and pressure
                        coo1 = Coordinates.fromGeodetic(rds(s).lat / 180 * pi, rds(s).lon / 180 * pi, [] ,rds(s).elevation);
                        launch_id = [1; (find(diff(datenum(rds(s).data_time)) > 0.4) + 1)];
                        launch_time = datenum(rds(s).data_time(launch_id));
                        pr1 = rds(s).pressure(launch_id);
                        
                        % Fill non valid values (< 6 sigma)
                        pr1((pr1 - median(pr1, 'omitnan') + 6 * std(pr1 - median(pr1, 'omitnan'))) < 0) = nan;
                        pr1 = interp1q(launch_time(~isnan(pr1)), pr1(~isnan(pr1)), launch_time);
                        
                        % Get GNSS station coordinates and pressure
                        coo2 = sta_list(id_rec(s)).getPos.getMedianPos;
                        [pr2, t2, h2] = sta_list(id_rec(s)).getPTH_mr;
                        time_pr2 = sta_list(id_rec(s)).getTime.getNominalTime.getMatlabTime;
                        
                        % find common epochs (around 30 min)
                        [~, idt_gnss, idt_raob] = intersect(round(time_pr2*48), round(launch_time*48));
                        
                        bias = pr2(idt_gnss) - pr1(idt_raob);
                        id_ko = pr1(idt_raob) == 1000; % it seems that all the values for a radiosonde less then 1000 (at ground) are saturated to 1000
                        if any(~id_ko)
                            bias = interp1q(launch_time(idt_raob(~id_ko)), bias(~id_ko), launch_time(idt_raob(id_ko)));
                            pr1(idt_raob(id_ko)) = pr2(idt_gnss(id_ko)) - bias;
                            % figure; plot(time_pr2(idt_gnss), pr2(idt_gnss), '*' ); hold on; plot(launch_time(idt_raob), pr1(idt_raob),'o'); setTimeTicks
                            
                            % This is my ZHD correction to move the closest
                            % station to the radiosonde data
                            time_corr = launch_time(idt_raob);
                            if Core.getState.isMet
                                zhd_corr = Atmosphere.getZenithDelayCorrection(coo1, pr1(idt_raob), coo2, pr2(idt_gnss));
                            else
                                zhd_corr = Atmosphere.getZenithDelayCorrection(coo1, pr1(idt_raob), coo2);
                            end
                        else
                            time_corr = time_pr2(idt_gnss);
                            % Use the data of the GNSS only -> radiosonde ppressure is not valid
                            zhd_corr = -Atmosphere.getZenithDelayCorrection(coo2, pr2(idt_gnss), coo1);
                        end
                        pause(0.1);
                        set(0, 'CurrentFigure', fh);
                        Core_Utils.plotSep(s_time.getMatlabTime, (s_ztd - interp1q(time_corr, zhd_corr, s_time.getMatlabTime)) * 1e2, '.-', 'LineWidth', 2);

                        [ztd_diff_sta, ztd_diff_sta_hcorr, id_mean] = deal(nan(time_rds.length, 1));
                        for e = 1 : time_rds.length
                            [t_min, id_min(e)] = min(abs(s_time - time_rds.getEpoch(e)));
                            if t_min < 3600
                                ztd_diff_sta(e) = ztd_rds(e) - 1e2.*s_ztd(id_min(e));                                
                                ztd_diff_sta_hcorr(e) = ztd_diff_sta(e) + 1e2 .* interp1q(time_corr, zhd_corr, time_pr2(id_min(e)));
                            end
                        end
                        m_diff_sta(s) = mean(ztd_diff_sta, 1, 'omitnan');
                        s_diff_sta(s) = std(ztd_diff_sta, 1, 'omitnan');
                        
                        m_diff_sta_hcorr(s) = mean(ztd_diff_sta_hcorr, 1, 'omitnan');
                        s_diff_sta_hcorr(s) = std(ztd_diff_sta_hcorr, 1, 'omitnan');
                        log.addMonoMessage(sprintf('%2d)  %6.2f cm    %6.2f cm      %4s  %9.1f   %9.1f   "%s"', s, m_diff_sta_hcorr(s), s_diff_sta_hcorr(s), sta_list(id_rec(s)).getMarkerName4Ch, round(d3d(s) / 1e3), dup(s), rds(s).getName()));
                        
                        pause(0.1);
                        set(0, 'CurrentFigure', fh);
                        plot(time_rds.getMatlabTime, ztd_rds, '.k', 'MarkerSize', 30, 'LineWidth', 2);
                        plot(time_rds.getMatlabTime, ztd_rds, '.w', 'MarkerSize', 40, 'LineWidth', 2);
                        plot(time_rds.getMatlabTime, ztd_rds, '.k', 'MarkerSize', 30, 'LineWidth', 2);
                        outm = {'ZTD GPS from interpolation', ...
                            sprintf('ZTD GPS of %s', sta_list(id_rec(s)).getMarkerName4Ch), ...
                            sprintf('ZTD GPS of %s @ RAOB height and pressure', sta_list(id_rec(s)).getMarkerName4Ch), ...
                            sprintf('ZTD RAOB @ %s', rds(s).getName())};
                        [~, icons] = legend(outm, 'location', 'southwest');
                        n_entry = numel(outm);
                        icons = icons(n_entry + 2 : 2 : end);
                        for i = 1 : numel(icons)
                            icons(i).MarkerSize = 18;
                        end
                        title(sprintf('Summary comparison @ %d Km, %.1f m up (mean, std)\n %s station ZHD corrected (%.2f cm, %.2f cm)\n %s station (%.2f cm, %.2f cm)\n Interpolated (%.2f cm, %.2f cm)\\fontsize{5} \n', ...
                            round(d3d(s) / 1e3), dup(s), ...
                            sta_list(id_rec(s)).getMarkerName4Ch, m_diff_sta_hcorr(s), s_diff_sta_hcorr(s), ...
                            sta_list(id_rec(s)).getMarkerName4Ch, m_diff_sta(s), s_diff_sta(s), ...
                            m_diff(s), s_diff(s)), 'FontName', 'Open Sans');
                        
                        xlim(s_time.getEpoch([1, s_time.length]).getMatlabTime);
                        ylabel(sprintf('ZTD [cm] comparison @ %d Km (%.1f m up)\\fontsize{5} \n', round(d3d(s) / 1e3), dup(s)));
                        setTimeTicks; grid minor;
                        drawnow;
                        ax = gca; ax.FontSize = 16;
                        Core_UI.beautifyFig(fh);
                        Core_UI.addExportMenu(fh);             
                        Core_UI.addBeautifyMenu(fh);             
                        fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                        
                        % % Histogram
                        % pause(0.1);
                        % fh = figure('Visible', 'off');
                        % Core_UI.beautifyFig(fh);
                        % fh.Name = sprintf('%03d: Rds Hist %s', fh.Number, rds(s).sta_num); fh.NumberTitle = 'off';                        
                        % fh_list = [fh_list; fh]; %#ok<AGROW>
                        % fig_name = sprintf('Raob_ZTD_%s_validation_hist', rds(s).sta_num);
                        % fh.UserData = struct('fig_name', fig_name);
                        % hist(ztd_diff_sta_hcorr,25);
                        % title('Histogram of ZTD RAOB - GNSS (ZHD corrected)');
                        % xlabel('dZTD [cm]');
                        % Core_UI.beautifyFig(fh);
                        % Core_UI.addExportMenu(fh);             
                        % Core_UI.addBeautifyMenu(fh);             
                        % fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;                        
                    end                     
                end
                
                % Plot map of all the radiosondes tests ----------------------------------------------
                % Retrieve DTM model
                
                if flag_show
                    log.addMarkedMessage('Preparing map, please wait...');
                    
                    % Radiometers points
                    data_mean = m_diff_sta_hcorr;
                    data_std = s_diff_sta_hcorr;
                    data_lat = rds.getLat();
                    data_lon = rds.getLon();

                    fh = sta_list.showMapDtm([], [], data_lat, data_lon);
                    fh.UserData.fig_name = 'Rds_Validation_Map';
                    fh.Name = sprintf('%03d: Raob ZTD Val.', fh.Number); fh.NumberTitle = 'off';
                    fh_list = [fh_list; fh];
                    gnss_list = findall(fh.Children(end).Children, 'UserData', 'GNSS_point');
                    for r = 1 : numel(gnss_list)
                        gnss_list(r).Color = Core_UI.LBLUE;
                    end
                    
                    n_col = round(max(abs(minMax(data_mean))*10));
                    caxis(n_col * [0 1] ./ 10); colormap(Cmap.get('linspaced', n_col));
                    drawnow
                    try
                        % It seems that the only way to have a colorbar correctly moving with the figure is to use the internal colorbar object
                        cb = colorbar;
                        %cb_m = m_contfbar(0.97, cb.Position(2) + [0 cb.Position(4)],[0 n_col/10], 0:0.1:(n_col/10),'edgecolor','none','endpiece','no', 'fontsize', 16);
                        %cb.Units = 'pixels';                        
                        %cb_m.Units = 'pixels';                        
                        %cb_m.Position(1) = cb.Position(1) + 10;
                        %ax_pos = ax.Position;
                        %delete(cb); % deleting the colorbar changes the size of the axes
                        %ax.Position = ax_pos;
                        %cb_m.Units = 'normalized';
                        xlabel(cb, 'cm','color','k');
                    catch
                        drawnow
                    end
                                        
                    [x, y] = m_ll2xy(data_lon, data_lat);
                    set(0, 'CurrentFigure', fh);
                    plot(x(:), y(:),'.k', 'MarkerSize', 5);
                    % Label BG (in background w.r.t. the point)
                    for r = 1 : numel(rds)
                        %name = rds(r).getName;
                        name = sprintf('%.1f, %.1f', data_mean(r), data_std(r));
                        text(x(r), y(r), ['     ' name ' '], ...
                            'FontWeight', 'bold', 'FontSize', 12, 'Color', [1 1 1], ...
                            'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                            'Margin', 2, 'LineWidth', 2, ...
                            'HorizontalAlignment','left');
                    end
                    
                    for r = 1 : numel(rds)
                        %name = rds(r).getName;
                        name = sprintf('%.1f, %.1f', data_mean(r), data_std(r));
                        t = text(x(r), y(r), ['     ' name ' '], ...
                            'FontWeight', 'bold', 'FontSize', 12, 'Color', [0 0 0], ...
                            'Margin', 2, 'LineWidth', 2, ...
                            'HorizontalAlignment','left');
                    end
                                        
                    set(0, 'CurrentFigure', fh);
                    %col_data = Cmap.getColor(round(data_mean * 10) + n_col, 2 * n_col, 'RdBu');
                    col_data = Cmap.getColor(max(1, round(abs(data_mean) * 10) + 1), n_col + 1, 'linspaced');
                    plot(x, y, 's', ...
                        'MarkerSize', 28, ...
                        'MarkerFaceColor', [0 0 0], ...
                        'MarkerEdgeColor', [0 0 0], ...
                        'UserData', 'RAOB_point');
                    for r = 1 : numel(rds)
                        plot(x(r), y(r), 's', ...
                            'MarkerSize', 24, ...
                            'MarkerFaceColor', col_data(r,:), ...
                            'UserData', 'RAOB_point');
                    end
                    
                    title(sprintf('Map of Radiosonde Validation\nUsing the closest GNSS Station (ZHD corrected)\nMean [cm] and StD [cm]\\fontsize{5} \n'), 'FontSize', 16);
                    Core_UI.beautifyFig(fh);
                    Core_UI.addExportMenu(fh);
                    Core_UI.addBeautifyMenu(fh);
                end
            end
        end
        
        function fh_list = showIgsValidation(sta_list)
            % Function to show the comparison between results stored in the
            % rec and official IGS solutions
            %
            % SYNTAX
            %   sta_list.showIgsPosValidation()
            
            log = Core.getLogger();
            state = Core.getCurrentSettings;
            fh_list = [];
            n_rec = length(sta_list);
            east_stat = nan(n_rec,2);
            north_stat = nan(n_rec,2);
            up_stat = nan(n_rec,2);
            ztd_stat = nan(n_rec,2);
            gn_stat = nan(n_rec,2);
            ge_stat = nan(n_rec,2);
            sta_names = {};
            id_ok = true(n_rec, 1);
            
            %[missing_days] = this.downloadIGSOfficialStation(sta_name, time, flag_download);
            
            for r = 1:n_rec
                flag_plot = true;
                if ~sta_list(r).out.isEmpty
                    xyz_igs = sta_list(r).out.getIgsXYZ;
                    %coo_gogps = Coordinates.fromXYZ(sta_list(r).out.xyz, sta_list(r).out.time_pos);
                    %coo_gogps.addOffset(-[sta_list(r).ant_delta_en sta_list(r).ant_delta_h]);       % Add antennna offset
                    %xyz_gogps = coo_gogps.getXYZ;
                    xyz_gogps = sta_list(r).out.getPos.getXYZ;
                    log.addMarkedMessage(sprintf('Analysing station "%s"', sta_list(r).getMarkerName4Ch));
                    if ~any(xyz_igs(:))
                        id_ok(r) = false;
                        log.addError(sprintf('"%s" does not appear in the IGS combination of AC global solutions', sta_list(r).getMarkerName4Ch));
                    else
                        coo_gogps = sta_list(r).getPos;
                        coo_igs = Coordinates.fromXYZ(xyz_igs, coo_gogps.time);
                        if flag_plot && coo_igs.time.length > 1
                            fh = Coordinates.showCompareENU([coo_igs; coo_gogps]);
                            set(0, 'CurrentFigure', fh);
                            ax = subplot(3,1,1);
                            ax.Title.String{1} = [sta_list(r).getMarkerName4Ch  ax.Title.String{1}(9:end)];
                            % If I have only one receiver use as name the name of the receiver
                            fig_name = sprintf('IGS_vs_goGPS_%s_%s', sta_list(r).getMarkerName4Ch, state.getSessionsStart.toString('yyyymmdd_HHMM'));                        
                            fh.UserData = struct('fig_name', fig_name);
                            Core_UI.addExportMenu(fh);
                            fh_list = [fh_list; fh];   
                            drawnow
                        end

                        xyz_diff = xyz_gogps - xyz_igs;
                        if ~any(xyz_diff(:))
                            log.addError(sprintf('"%s" does not have valid positions for goGPS or IGS', sta_list(r).getMarkerName4Ch));
                            id_ok(r) = false;
                        else
                            enu_diff = Coordinates.cart2local(sta_list(r).out.getMedianPosXYZ, xyz_diff);
                            sensor = enu_diff - repmat(median(enu_diff,'omitnan'),size(enu_diff,1),1);
                            out_idx = sum((abs(sensor) > 0.05),2) > 0;
                            enu_diff(out_idx,:) =[];
                            east_stat(r,:) = [mean(enu_diff(:,1),'omitnan'),perc(abs(enu_diff(:,1) - mean(enu_diff(:,1),'omitnan')),0.95)]*1e3;
                            north_stat(r,:) = [mean(enu_diff(:,2),'omitnan'),perc(abs(enu_diff(:,2) - mean(enu_diff(:,2),'omitnan')),0.95)]*1e3;
                            up_stat(r,:) = [mean(enu_diff(:,3),'omitnan'),perc(abs(enu_diff(:,3) - mean(enu_diff(:,3),'omitnan')),0.95)]*1e3;
                            [ztd_diff, gn_diff, ge_diff] =  sta_list(r).out.getIgsTropo('interp_value');
                            ztd_diff = ztd_diff - sta_list(r).out.ztd;
                            gn_diff = gn_diff - sta_list(r).out.tgn;
                            ge_diff = ge_diff - sta_list(r).out.tge;
                            out_idx = abs(ztd_diff) >0.05 | abs(gn_diff) >0.01 | abs(ge_diff) > 0.01;
                            ztd_diff(out_idx) = [];
                            gn_diff(out_idx) = [];
                            ge_diff(out_idx) = [];
                            ztd_stat(r,:) = [mean(ztd_diff,'omitnan'),perc(abs(ztd_diff - mean(ztd_diff,'omitnan')),0.95)]*1e3;
                            gn_stat(r,:) = [mean(gn_diff,'omitnan'),perc(abs(gn_diff - mean(gn_diff,'omitnan')),0.95)]*1e3;
                            ge_stat(r,:) = [mean(ge_diff,'omitnan'),perc(abs(ge_diff - mean(ge_diff,'omitnan')),0.95)]*1e3;
                        end
                    end
                end
                sta_names{end+1} = lower(sta_list(r).getMarkerName4Ch);
            end
            if any(id_ok)             
                % sort by bet on the east axis
                %[~,idx] = sort(abs(east_stat(:,1)));
                [~, idx] = sort(east_stat(:,1).^2 + north_stat(:,1).^2 + 0*up_stat(:,1).^2);
                idx = idx(id_ok(idx));
                east_stat = east_stat(idx,:);
                north_stat = north_stat(idx,:);
                up_stat = up_stat(idx,:);
                ztd_stat = ztd_stat(idx,:);
                gn_stat = gn_stat(idx,:);
                ge_stat = ge_stat(idx,:);
                sta_names = sta_names(idx);
                
                sta_list = sta_list(idx);
                
                n_rec = length(sta_list);
                
                fh = figure('Visible', 'off');
                fh.Name = sprintf('%03d: IGS Validation', fh.Number); fh.NumberTitle = 'off';
                
                fh_list = [fh_list; fh];
                if numel(sta_list) == 1
                    % If I have only one receiver use as name the name of the receiver
                    fig_name = sprintf('IGS_Comparison_%s_%s', sta_list.getMarkerName4Ch, state.getSessionsStart.toString('yyyymmdd_HHMM'));
                else
                    % If I have more than one receiver use as name the name of the project
                    fig_name = sprintf('IGS_Comparison_%s', strrep(Core.getState.getPrjName,' ', '_'));
                end
                fh.UserData = struct('fig_name', fig_name);
                
                subplot(3,2,1)
                errorbar(1:n_rec,east_stat(:,1),east_stat(:,2),'.','MarkerSize',25,'LineWidth',2,'Color',Core_UI.getColor(1,6))
                ylabel('[mm]')
                title('Position East')
                ax = gca;
                ax.YGrid = 'on';
                ax.GridLineStyle = '-';
                set(gca, 'YTick', [-30 -10 0 10 30])
                ylim([-30 30])
                xlim([0.5, n_rec + 0.5]);
                set(gca, 'XTickLabels', {})
                set(gca,'fontweight','bold','fontsize',16)
                setAllLinesWidth(1.3)
                
                subplot(3,2,3)
                errorbar(1:n_rec,north_stat(:,1),north_stat(:,2),'.','MarkerSize',25,'LineWidth',2,'Color',Core_UI.getColor(2,6))
                ax = gca;
                ax.YGrid = 'on';
                ax.GridLineStyle = '-';
                set(gca, 'YTick', [-30 -10 0 10 30])
                ylim([-30 30])
                xlim([0.5, n_rec + 0.5]);
                set(gca, 'XTickLabels', {})
                ylabel('[mm]')
                set(gca,'fontweight','bold','fontsize',16)
                title('Position North')
                
                subplot(3,2,5)
                errorbar(1:n_rec,up_stat(:,1),up_stat(:,2),'.','MarkerSize',25,'LineWidth',2,'Color',Core_UI.getColor(3,6))
                ax = gca;
                ax.YGrid = 'on';
                ax.GridLineStyle = '-';
                ylim([-50 50]);
                xlim([0.5, n_rec + 0.5]);
                set(gca, 'YTick', [-50 -20 0 20 50])
                %set(gca, 'XTick', [1:28])
                set(gca, 'XTickLabels', {})
                ylabel('[mm]')
                title('Position Up')
                set(gca, 'XTickLabelRotation', 45)
                set(gca,'fontweight','bold','fontsize',16)
                
                set(gca, 'XTick', [1:n_rec]);
                set(gca, 'XTickLabels', upper(sta_names))
                set(gca, 'XTickLabelRotation', 45)
                set(gca,'fontweight','bold','fontsize',16)
                
                subplot(3,2,2)
                errorbar(1:n_rec,ztd_stat(:,1),ztd_stat(:,2),'.','MarkerSize',25,'LineWidth',2,'Color',Core_UI.getColor(4,6))
                ylabel('[mm]')
                xlim([0.5, n_rec + 0.5]);
                set(gca, 'XTickLabels', {})
                ylim([-25 25])
                ax = gca;
                
                ax.YGrid = 'on';
                ax.GridLineStyle = '-';
                set(gca, 'YTick', [-25 -10 0 10 25])
                set(gca,'fontweight','bold','fontsize',16)
                title('ZTD')
                
                subplot(3,2,4)
                errorbar(1:n_rec,gn_stat(:,1),gn_stat(:,2),'.','MarkerSize',25,'LineWidth',2,'Color',Core_UI.getColor(5,6))
                ylabel('[mm]')
                set(gca, 'XTickLabels', {})
                ylim([-4 4])
                xlim([0.5, n_rec + 0.5]);
                ax = gca;
                
                ax.YGrid = 'on';
                ax.GridLineStyle = '-';
                set(gca, 'YTick', [-4 -1 0 1 4])
                set(gca,'fontweight','bold','fontsize',16)
                title('North gradient')
                
                subplot(3,2,6)
                errorbar(1:n_rec,ge_stat(:,1),ge_stat(:,2),'.','MarkerSize',25,'LineWidth',2,'Color',Core_UI.getColor(6,6))
                ylabel('[mm]')
                ylim([-4 4])
                xlim([0.5, n_rec + 0.5]);
                ax = gca;
                
                ax.YGrid = 'on';
                ax.GridLineStyle = '-';
                set(gca, 'YTick', [-4 -1 0 1 4])
                set(gca, 'XTickLabels', {})
                title('East gradient');
                
                set(gca, 'XTick', [1:n_rec]);
                set(gca, 'XTickLabels', upper(sta_names))
                set(gca, 'XTickLabelRotation', 45)
                set(gca,'fontweight','bold','fontsize',16)
                
                Core_UI.addExportMenu(fh);
                Core_UI.addBeautifyMenu(fh);
                Core_UI.beautifyFig(fh, 'dark');
                fh.Visible = iif(Core_UI.isHideFig, 'off', 'on');
            end
        end
        
        function [fh_list, m_diff, s_diff] = showIgsTropoValidation(sta_list, igs_list)
            % Compute and show comparison with igs station
            % given region list, and station id (as cell arrays)
            %
            % INPUT
            %   sta_list        list of gnss receivers
            %   rds_list        cell array of string containing the ig station ID
            %
            % OUTPUT
            %   m_diff          mean of the ZTD differences
            %   s_diff          std of the ZTD differences
            %
            % SYNTAX
            %  [m_diff, s_diff] = sta_list.showIgsTropoValidation(station_list);
            %
            % EXAMPLE
            %  % testing geonet full network
            %  [m_diff, s_diff] = sta_list.showIgsTropoValidation(Radiosonde.JAPAN_STATION);
            %
            % SEE ALSO
            %   Radiosonde GNSS_Station.getIgsTropoValidation
            if nargin == 1
                igs_list = [];
            end
            sta_list = sta_list(~sta_list.isEmptyOut_mr);
            [m_diff, s_diff, fh_list] = sta_list.getIgsTropoValidation(igs_list, true);
        end
        
        function [m_diff, s_diff, fh_list] = getIgsTropoValidation(sta_list, igs_list, flag_show, flag_download)
            % Compute a comparison with IGS ststion
            % given region list, and station id (as cell arrays)
            %
            % INPUT
            %   sta_list        list of gnss receivers
            %   igs_list        cell array of string containing the IGS station ID 
            %
            % OUTPUT
            %   m_diff          mean of the ZTD differences
            %   s_diff          std of the ZTD differences
            %
            % SYNTAX
            %  [m_diff, s_diff] = sta_list.getIgsTropoValidation(station_list);
            %
            % EXAMPLE
            %  % testing geonet full network
            %  [m_diff, s_diff] = sta_list.getIgsTropoValidation('usno');
            %
            
            if nargin < 4 || isempty(flag_download)
                flag_download = true;
            end
            
            sta_list = sta_list(~sta_list.isEmptyOut_mr);
            log = Core.getLogger();
            state = Core.getCurrentSettings;
            [m_diff, s_diff, fh_list] = deal([]);
            if nargin < 3
                 flag_show = false;
            end
            
            if nargin < 2 || isempty(igs_list)
                igs_list = sta_list.getMarkerName4Ch();
            end
            
            log.addMarkedMessage('Retrieving data, please wait...');
            % Get time limits
            p_time = sta_list.getSyncTimeExpanded(sta_list);
            start_time = GPS_Time(floor(p_time.first.getMatlabTime * 2)/2);
            stop_time = GPS_Time(ceil(p_time.last.getMatlabTime * 2)/2);

            % Download IGS Stations
            tsc = Tropo_Sinex_Compare();
            tsc.addIGSOfficialStation(igs_list, p_time, flag_download);
            
            if isempty(tsc.results) || isempty(fieldnames(tsc.results))
                log.addWarning('No IGS valid station found for comparison');
            else
                % For now disable interpolation, check only IGS stations
                % log.addMarkedMessage('Get GNSS interpolated ZTD @ IGS locations');
                % [ztd, ztd_height_correction, time] = sta_list.getTropoInterp('ZTD', tsc.getLat(), tsc.getLon(), tsc.getHeightOrtho(), 300);
                
                % Get closer GNSS stations
                [id_rec, d3d, dup] = sta_list.getCloserRec(tsc.getLat(), tsc.getLon(), tsc.getHeightOrtho());
                gnss_list = sta_list(id_rec);
                
                % Compute values
                log.addMonoMessage(sprintf('---------------------------------------------------------------------------------------'));
                [m_diff, s_diff] = deal(nan(tsc.getNumberSinex(), 1));
                log.addMonoMessage(sprintf(' ZTD IGS Validation\n---------------------------------------------------------------------------------------'));
                log.addMonoMessage(sprintf('                                Marker              Elevation     '));
                log.addMonoMessage(sprintf('       Mean          Std         GNSS    Dist [km]  diff. [m]  IGS Station'));
                log.addMonoMessage(sprintf('---------------------------------------------------------------------------------------'));
                
                % For each IGS station
                for s = 1 : tsc.getNumberSinex()
                    % radiosondes
                    [ztd_igs, time_igs] = tsc.getZtdSinex(s);
                    [gn_igs, ge_igs] = tsc.getGradientsSinex(s);
                    ztd_igs = ztd_igs * 100;
                    % gnss
                    [ztd_gogps, time_gogps] = gnss_list(s).getZtd;
                    gn_gogps = gnss_list(s).getTropoPar('gn');
                    ge_gogps = gnss_list(s).getTropoPar('ge');
                    
                    [~, id_igs, id_gogps] = intersect(round(time_igs.getMatlabTime * (86400 / 300)), round(time_gogps.getMatlabTime * (86400 / 300)));
                    
                    ztd_diff = ztd_igs(id_igs) - ztd_gogps(id_gogps)*1e2;
                    gn_diff = gn_igs(id_igs)*1e3 - gn_gogps(id_gogps)*1e3;
                    ge_diff = ge_igs(id_igs)*1e3 - ge_gogps(id_gogps)*1e3;
                    
                    m_diff(s) = mean(ztd_diff, 1, 'omitnan');
                    s_diff(s) = std(ztd_diff, 1, 'omitnan');
                    log.addMonoMessage(sprintf('%2d)  %6.2f cm    %6.2f cm      %4s  %9.1f   %9.1f   "%s"\n', s, m_diff(s), s_diff(s), gnss_list(s).getMarkerName4Ch, round(d3d(s) / 1e3), dup(s), tsc.getName(s)));
                    
                    if flag_show
                        fh = figure('Visible', 'on'); Core_UI.beautifyFig(fh);
                        fh_list = [fh_list; fh]; %#ok<AGROW>
                        subplot(2,2,1);
                        fh.Name = sprintf('%03d: dZTD %s', fh.Number, gnss_list(s).getMarkerName4Ch); fh.NumberTitle = 'off';
                        fh_list = [fh_list; fh]; %#ok<AGROW>
                        fig_name = sprintf('IGS_dZTD_%s_Validation', gnss_list(s).getMarkerName4Ch);
                        fh.UserData = struct('fig_name', fig_name);
                        Core_Utils.plotSep(time_igs.getEpoch(id_igs).getMatlabTime, ztd_diff, '.-');
                        xlim(time_igs.getEpoch(id_igs([1 end])).getMatlabTime);
                        ylabel('dZTD [cm]')
                        title(sprintf('%s IGS vs. goGPS [ cm ]\\fontsize{5} \n', gnss_list(s).getMarkerName4Ch));
                        setTimeTicks(3);
                        
                        subplot(2,2,2);
                        hist(ztd_diff,100);
                        title(sprintf('%s IGS vs. goGPS [ cm ]\ndZTD  Mean %.2f cm, std %.2fcm\\fontsize{5} \n', gnss_list(s).getMarkerName4Ch, m_diff(s), s_diff(s)));
                        
                        subplot(2,2,3);
                        hist(gn_diff,100);
                        title(sprintf('%s IGS vs. goGPS [ mm ]\ndGN Mean %.2f mm, std %.2fmm\\fontsize{5} \n', gnss_list(s).getMarkerName4Ch, mean(gn_diff, 1, 'omitnan'), mean(gn_diff, 1, 'omitnan')));
                        
                        subplot(2,2,4);
                        hist(ge_diff,100);
                        title(sprintf('%s IGS vs. goGPS [ mm ]\ndGE Mean %.2f mm, std %.2fmm\\fontsize{5} \n', gnss_list(s).getMarkerName4Ch, mean(ge_diff, 1, 'omitnan'), mean(ge_diff, 1, 'omitnan')));
                        
                        Core_UI.addExportMenu(fh);
                        Core_UI.addBeautifyMenu(fh);
                        Core_UI.beautifyFig(fh);
                    end
                    
                    coo_gogps = gnss_list(s).getPos;
                    if isfield(tsc.results.r2, (['r' gnss_list(s).getMarkerName4Ch]))
                        coo_igs = Coordinates.fromXYZ(tsc.results.r2.(['r' gnss_list(s).getMarkerName4Ch]).xyz, tsc.results.r2.(['r' gnss_list(s).getMarkerName4Ch]).coord_time);
                        
                        fh = Coordinates.showCompareENU([coo_igs, coo_gogps]);
                        set(0, 'CurrentFigure', fh);
                        ax = subplot(3,1,1);
                        ax.Title.String{1} = [gnss_list(s).getMarkerName4Ch  ax.Title.String{1}(9:end)];
                        fig_name = sprintf('IGS_TropoCoo_vs_goGPS_%s_%s', gnss_list(s).getMarkerName4Ch, state.getSessionsStart.toString('yyyymmdd_HHMM'));
                        fh.UserData = struct('fig_name', fig_name);
                        Core_UI.addExportMenu(fh);
                        fh_list = [fh_list; fh]; %#ok<AGROW>
                    end
                end
                log.addMonoMessage(sprintf('---------------------------------------------------------------------------------------\n'));
                
                if flag_show
                    % Plot comparisons
                    for s = 1 : tsc.getNumberSinex()
                        %%
                        fh = figure('Visible', 'off');
                        fh.Name = sprintf('%03d: IGS Tropo %s', fh.Number, gnss_list(s).getMarkerName4Ch); fh.NumberTitle = 'off';
                        
                        fh_list = [fh_list; fh]; %#ok<AGROW>
                        fig_name = sprintf('IGS_Tropo_vs_goGPS_%s_%s', gnss_list(s).getMarkerName4Ch, state.getSessionsStart.toString('yyyymmdd_HHMM'));
                        fh.UserData = struct('fig_name', fig_name);
                        
                        % interpolated ZTD
                        %Core_Utils.plotSep(time.getMatlabTime, ztd(:,s) + ztd_height_correction(s), '.-', 'LineWidth', 2, 'MarkerSize', 5);
                        %hold on;
                        
                        % closer ZTD
                        [s_ztd, s_time] = gnss_list(s).getZtd;
                        Core_Utils.plotSep(s_time.getMatlabTime, s_ztd * 1e2, '.-', 'LineWidth', 2, 'MarkerSize', 5); hold on;
                        
                        % IGS values
                        [ztd_igs, time_igs] = tsc.getZtdSinex(s);
                        Core_Utils.plotSep(time_igs.getMatlabTime, ztd_igs*100, '.-k', 'MarkerSize', 5, 'LineWidth', 2);
                        %outm = {'ZTD GPS from interpolation', sprintf('ZTD GPS of %s', sta_list(id_rec(s)).getMarkerName4Ch), ...
                        %    sprintf('IGS ZTD @ %s', tsc.getName(s))};
                        outm = {sprintf('ZTD GPS of %s', gnss_list(s).getMarkerName4Ch), ...
                            sprintf('IGS ZTD @ %s', tsc.getName(s))};
                        [~, icons] = legend(outm, 'location', 'northwest');
                        n_entry = numel(outm);
                        icons = icons(n_entry + 2 : 2 : end);
                        for i = 1 : numel(icons)
                            icons(i).MarkerSize = 18;
                        end
                        if (d3d(s) < 2.5)
                            title(sprintf('%s IGS comparison [cm]\\fontsize{5} \n', gnss_list(s).getMarkerName4Ch));
                        else
                            title(sprintf('%s IGS comparison @ %d Km (%.1f m up)\\fontsize{5} \n', gnss_list(s).getMarkerName4Ch, round(d3d(s) / 1e3), dup(s)));
                        end
                        ylabel('ZTD [cm]');
                        xlim(minMax(time_igs.getMatlabTime));
                        setTimeTicks; grid minor;
                        drawnow;
                        ax = gca; ax.FontSize = 16;
                        Core_UI.beautifyFig(gcf);
                        Core_UI.addExportMenu(gcf);
                        Core_UI.addBeautifyMenu(gcf);
                        fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                    end
                end
                
                % Plot map of all the radiosondes tests ----------------------------------------------
                % Retrieve DTM model
                
                if flag_show
                    Core.getLogger.addMarkedMessage('Preparing map, please wait...');
                    % set map limits
                    lat = tsc.getLat;
                    lon = tsc.getLon;
                    % set map limits
                    
                    if numel(lon) == 1
                        lon_lim = minMax(lon) + [-0.05 0.05];
                        lat_lim = minMax(lat) + [-0.05 0.05];
                    else
                        lon_lim = minMax(lon); lon_lim = min(179.999, max(-179.999, lon_lim + [-1 1] * diff(lon_lim) / 6));
                        lat_lim = minMax(lat); lat_lim = min(89.999, max(-89.999,lat_lim + [-1 1] * diff(lat_lim) / 6));
                    end
                    lon_lim = lon_lim + [-1 1] * max(0, (0.5 - diff(minMax(lon_lim))) / 2);
                    lat_lim = lat_lim + [-1 1] * max(0, (0.5 - diff(minMax(lat_lim))) / 2);
                    nwse = [lat_lim(2), lon_lim(1), lat_lim(1), lon_lim(2)];
                    clon =  min(180, max(-180, nwse([2 4]) + [-0.001 0.001]));
                    clat = max(-90, min(90, nwse([3 1]) + [-0.001 0.001]));
            
                    fh = figure('Visible', 'off');
                    fh.Color = [1 1 1];
                    
                    fh_list = [fh_list; fh];
                    fig_name = 'IGS_Validation_map';
                    fh.UserData = struct('fig_name', fig_name);
                    fh.Name = sprintf('%03d: IGS ZTD Val.', fh.Number); fh.NumberTitle = 'off';
                    
                    m_proj('equidistant','lon',clon,'lat',clat);   % Projection
                    % m_proj('utm', 'lon',clon,'lat',clat);   % Projection
                    axes
                    cmap = flipud(gray(1000));
                    colormap(cmap(150: end, :));
                    
                    % retrieve external DTM
                    try
                        [dtm, lat_dtm, lon_dtm] = Core.getRefDTM(nwse, 'ortho', 'low');
                        dtm = flipud(dtm);
                        % comment the following line to have bathimetry
                        dtm(dtm < 0) = nan; % - 1/3 * max(dtm(:));
                        
                        % uncomment the following line to have colors
                        % colormap(Cmap.adaptiveTerrain(minMax(dtm(:))));
                        drawnow;
                        
                        [shaded_dtm] = m_shadedrelief(lon_dtm, lat_dtm, dtm, 'nan', [0.98, 0.98, 1], 'gra', 30);
                        
                        %h_dtm = m_pcolor(lon_dtm, lat_dtm, dtm);
                        %h_dtm.CData = shaded_dtm;
                        m_image(lon_dtm, lat_dtm, shaded_dtm);
                    catch
                        % use ETOPO1 instead
                        m_etopo2('shadedrelief','gradient', 3);
                    end
                    
                    % read shapefile
                    shape = 'none';
                    if (~strcmp(shape,'none'))
                        if (~strcmp(shape,'coast')) && (~strcmp(shape,'fill'))
                            if (strcmp(shape,'10m'))
                                M = m_shaperead('countries_10m');
                            elseif (strcmp(shape,'30m'))
                                M = m_shaperead('countries_30m');
                            else
                                M = m_shaperead('countries_50m');
                            end
                            [x_min, y_min] = m_ll2xy(min(lon_lim), min(lat_lim));
                            [x_max, y_max] = m_ll2xy(max(lon_lim), max(lat_lim));
                            for k = 1 : length(M.ncst)
                                lam_c = M.ncst{k}(:,1);
                                ids = lam_c <  min(lon);
                                lam_c(ids) = lam_c(ids) + 360;
                                phi_c = M.ncst{k}(:,2);
                                [x, y] = m_ll2xy(lam_c, phi_c);
                                if sum(~isnan(x))>1
                                    x(find(abs(diff(x)) >= abs(x_max - x_min) * 0.90) + 1) = nan; % Remove lines that occupy more than th 90% of the plot
                                    line(x,y,'color', [0.3 0.3 0.3]);
                                end
                            end
                        else
                            if (strcmp(shape,'coast'))
                                m_coast('line','color', lineCol);
                            else
                                m_coast('patch',lineCol);
                            end
                        end
                    end
                    
                    hold on;
                    m_grid('box','fancy','tickdir','in', 'fontsize', 16);
                    
                    drawnow
                    m_ruler([.7 1], -0.05, 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
                    
                    % IGS points
                    data_mean = m_diff;
                    data_std = s_diff;
                    data_lat = tsc.getLat();
                    data_lon = tsc.getLon();
                    
                    [x, y] = m_ll2xy(data_lon, data_lat);
                    
                    plot(x(:), y(:),'.k', 'MarkerSize', 5);
                    
                    n_col = round(max(abs(minMax(data_mean))*10)) + 1;
                    caxis(n_col * [0 1] ./ 10); colormap(Cmap.get('linspaced', n_col));
                    drawnow
                    try
                        % It seems that the only way to have a colorbar correctly moving with the figure is to use the internal colorbar object
                        cb = colorbar;
                        %cb_m = m_contfbar(0.97, cb.Position(2) + [0 cb.Position(4)],[0 n_col/10], 0:0.1:(n_col/10),'edgecolor','none','endpiece','no', 'fontsize', 16);
                        %cb.Units = 'pixels';
                        %cb_m.Units = 'pixels';
                        %cb_m.Position(1) = cb.Position(1) + 10;
                        %ax_pos = ax.Position;
                        %delete(cb); % deleting the colorbar changes the size of the axes
                        %ax.Position = ax_pos;
                        %cb_m.Units = 'normalized';
                        xlabel(cb, 'cm','color','k');
                        drawnow
                    catch
                        drawnow
                    end
                    
                    % Label BG (in background w.r.t. the point)
                    for r = 1 : numel(gnss_list)
                        name = sprintf('%.2f, %.2f', data_mean(r), data_std(r));
                        text(x(r), y(r), ['     ' name ' '], ...
                            'FontWeight', 'bold', 'FontSize', 12, 'Color', [1 1 1], ...
                            'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                            'Margin', 2, 'LineWidth', 2, ...
                            'HorizontalAlignment','left');
                    end
                    
                    for r = 1 : tsc.getNumberSinex
                        name = sprintf('%.2f, %.2f', data_mean(r), data_std(r));
                        text(x(r), y(r), ['     ' name ' '], ...
                            'FontWeight', 'bold', 'FontSize', 12, 'Color', [0 0 0], ...
                            'Margin', 2, 'LineWidth', 2, ...
                            'HorizontalAlignment','left');
                    end
                    
                    %col_data = Cmap.getColor(round(data_mean * 10) + n_col, 2 * n_col, 'RdBu');
                    col_data = Cmap.getColor(floor(abs(data_mean) * 10) + 1, n_col, 'linspaced');
                    plot(x(:), y(:), '.', 'MarkerSize', 100, 'Color', [0 0 0]);
                    for r = 1 : numel(gnss_list)
                        plot(x(r), y(r), '.', 'MarkerSize', 92, 'Color', col_data(r,:));
                    end
                    
                    title(sprintf('Map of IGS Tropo Validation\nMean [cm] and StD [cm]\\fontsize{5} \n'), 'FontSize', 16);
                    
                    Core_UI.addExportMenu(fh);
                    Core_UI.addBeautifyMenu(fh);
                    Core_UI.beautifyFig(fh);
                    fh.Visible = iif(Core_UI.isHideFig, 'off', 'on');
                    Core.getLogger.addStatusOk('The map is ready ^_^');
                end
            end
        end        
    end
          
end
