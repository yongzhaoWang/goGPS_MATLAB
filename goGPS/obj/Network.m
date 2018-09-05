%   CLASS Network
% =========================================================================
%
% DESCRIPTION
%   Class to manage network processing of receivers
%
% EXAMPLE
%   Net = Network();
%
% SEE ALSO
% FOR A LIST OF CONSTANTs and METHODS use doc Network

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 4 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------
classdef Network < handle
    properties
        rec_list
        state
        
        common_time      % gps_time
        rec_time_indexes % indexes
        coo              % [n_coo x n_rec] receiver coordinates
        clock            % [n_epoch x n_rec] reciever clock
        ztd              % [n_epoch x n_rec] reciever ZTD
        ztd_gn           % [n_epoch x n_rec] reciever ZTD gradients north
        ztd_ge           % [n_epoch x n_rec] reciever ZTD gradients east
        amb              % {n_rec} recievers ambiguity
        log
    end
    methods
        function this = Network(rec_list)
            this.rec_list = rec_list;
            this.state = Global_Configuration.getCurrentSettings;
            this.log = Logger.getInstance();
        end
        
        function adjust(this, idx_ref)
            %  adjust the gnss network
            %
            % SYNATAX;
            %    this. adjustNetwork(idx_ref)
            % INPUT:
            %     idx_ref : [1,n_rec] boolean, receivers to be choosen as reference, their value mean will be set to zero
            
            % set up the the network adjustment
            if nargin < 2 || any(isnan(idx_ref))
                idx_ref = 1 : numel(this);
            end
            is_empty_recs = this.rec_list.isEmpty_mr;
            if sum(~is_empty_recs) > 1
                e = find(is_empty_recs);
                this.rec_list(e) = [];
                idx_ref(idx_ref == e) = [];
                ls = Least_Squares_Manipulator(this.rec_list(1).cc);
                [this.common_time, this.rec_time_indexes]  = ls.setUpNetworkAdj(this.rec_list);
                n_time = this.common_time.length;
                if this.state.flag_tropo
                    ls.setTimeRegularization(ls.PAR_TROPO, (this.state.std_tropo)^2 / 3600 * ls.rate );
                end
                if this.state.flag_tropo_gradient
                    ls.setTimeRegularization(ls.PAR_TROPO_N, (this.state.std_tropo_gradient)^2 / 3600 * ls.rate );
                    ls.setTimeRegularization(ls.PAR_TROPO_E, (this.state.std_tropo_gradient)^2 / 3600 * ls.rate );
                end
                [x, res] = ls.solve;
                s0 = mean(abs(res(res~=0)));
                this.log.addMessage(this.log.indent(sprintf('Network solution computed,  s0 = %.4f',s0)));
                if s0 < 0.01
                    % intilaize array for results
                    n_rec = length(this.rec_list);
                    this.clock = zeros(n_time, n_rec);
                    this.coo = nan(n_rec, 3);
                    this.ztd = nan(n_time, n_rec);
                    this.ztd_gn = nan(n_time, n_rec);
                    this.ztd_ge = nan(n_time, n_rec);
                    
                    % --- fill the correction values in the network
                    for i = 1 : n_rec
                        % if all value in the receiver are set to nan initilaize them to zero
                        if sum(isnan(this.rec_list(i).work.ztd)) == length(this.rec_list(i).work.ztd)
                            this.rec_list(i).work.ztd(:) = 0;
                            this.rec_list(i).work.tge(:) = 0;
                            this.rec_list(i).work.tgn(:) = 0;
                        end
                        % for all paramter take the apriori in the receiver and sum the netwrok estimated correction
                        idx_rec = x(:,3) == i;
                        coo = [x(x(:,2) == 1 & idx_rec,1) x(x(:,2) == 2 & idx_rec,1) x(x(:,2) == 3 & idx_rec,1)];
                        this.coo(i,:) = coo;
                        
                        clk = x(x(:,2) == ls.PAR_REC_CLK & idx_rec,1);
                        this.clock(~isnan(this.rec_time_indexes(:,i)),i) = clk;
                        
                        if this.state.flag_tropo
                            ztd = x(x(:,2) == ls.PAR_TROPO & idx_rec,1);
                            this.ztd(~isnan(this.rec_time_indexes(:,i)),i) = ztd;
                        end
                        
                        if this.state.flag_tropo_gradient
                            gn = x(x(:,2) == ls.PAR_REC_CLK & idx_rec,1);
                            this.ztd_gn(~isnan(this.rec_time_indexes(:,i)),i) = gn;
                            
                            ge = x(x(:,2) == ls.PAR_REC_CLK & idx_rec,1);
                            this.ztd_ge(~isnan(this.rec_time_indexes(:,i)),i) = ge;
                        end
                    end
                    
                    % ALL OF THIS MAKES NO SENSE TO ME (Andrea). Now it should ;) (Giulio)
                    %--- transform the result in the desired free network
                    
                    if ~isnan(idx_ref(1)) && ~(length(idx_ref) == 1 && idx_ref == 1)
                        
                        S = zeros(n_rec);
                        S(:, idx_ref) = - 1 / numel(idx_ref);
                        S = S + eye(n_rec);  % < - this should be an S trasform but i am not sure
                        % it is the paramter itself  the mean of the reference paramter
                        % it is in matrix form so it can be used in the future for variance covariance matrix of the coordinates
                        
                        % Applying the S transform I obtain the corrections with respect to the reference
                        this.coo(:,1) = S * this.coo(:,1);
                        this.coo(:,2) = S * this.coo(:,2);
                        this.coo(:,3) = S * this.coo(:,3);
                        
                        % apply the S transform to the epochwise parameters
                        for i = 1 : n_time
                            id_present = ~isnan(this.clock(i,:));
                            idx_ref_t = intersect(idx_ref, find(id_present));
                            if isempty(idx_ref_t)
                                S = nan;
                            else
                                n_rec_t = sum(id_present);
                                S = zeros(n_rec_t);
                                S(:,idx_ref) = - 1 / numel(idx_ref_t);
                                S = S + eye(n_rec_t);
                            end
                            % clock
                            this.clock(i,:) = (S*this.clock(i,:)')';
                            % ztd
                            if this.state.flag_tropo
                                this.ztd(i,:) = (S*this.ztd(i,:)')';
                            end
                            % gradients
                            if this.state.flag_tropo_gradient
                                this.ztd_gn(i,:) = (S*this.ztd_gn(i,:)')';
                                this.ztd_ge(i,:) = (S*this.ztd_ge(i,:)')';
                            end
                        end
                    end
                    
                    
                    % --- add the apriori values
                    for i = 1 : n_rec
                        % if all value in the receiver are set to nan initilaize them to zero
                        if sum(isnan(this.rec_list(i).work.ztd)) == length(this.rec_list(i).work.ztd)
                            this.rec_list(i).work.ztd(:) = 0;
                            this.rec_list(i).work.tge(:) = 0;
                            this.rec_list(i).work.tgn(:) = 0;
                        end
                        % for all paramter take the apriori in the receiver and sum the netwrok estimated correction
                        this.coo(i,:) = this.coo(i,:) + this.rec_list(i).work.xyz;
                        %
                        [idx_is, idx_pos] = ismembertol(this.rec_list(i).work.getTime.getGpsTime(), this.common_time.getGpsTime, 0.002, 'DataScale', 1);
                        idx_pos = idx_pos(idx_pos > 0);
                        clk_rec = this.rec_list(i).work.getDt();
                        this.clock(idx_pos,i) = this.clock(idx_pos,i) + clk_rec(idx_is);
                        
                        if this.state.flag_tropo
                            ztd_rec = this.rec_list(i).work.getZtd();
                            ztd_rec_apr = this.rec_list(i).work.getZwd() + this.rec_list(i).work.getAprZhd();
                            ztd_rec(ztd_rec == 0) = ztd_rec_apr(ztd_rec == 0);
                            this.ztd(idx_pos,i) = this.ztd(idx_pos,i) + ztd_rec(idx_is);
                        end
                        
                        if this.state.flag_tropo_gradient
                            [gn_rec, ge_rec] = this.rec_list(i).work.getGradient();
                            
                            this.ztd_gn(idx_pos,i) = this.ztd_gn(idx_pos,i) + gn_rec(idx_is);
                            
                            this.ztd_ge(idx_pos,i) = this.ztd_ge(idx_pos,i) + ge_rec(idx_is);
                        end
                    end
                    
                    % --- push back the results in the receivers
                    for i = 1 : n_rec
                        this.rec_list(i).work.xyz = this.coo(i,:);
                        idx_res_av = ~isnan(this.clock(:, i));
                        [idx_is, idx_pos] = ismembertol(this.common_time.getEpoch(idx_res_av).getGpsTime(), this.rec_list(i).work.time.getGpsTime, 0.002, 'DataScale', 1);
                        idx_pos = idx_pos(idx_pos > 0);
                        clk = this.clock(idx_res_av, i);
                        this.rec_list(i).work.dt(idx_pos) = clk(idx_is);
                        if this.state.flag_tropo
                            ztd = this.ztd(idx_res_av, i);
                            this.rec_list(i).work.ztd(idx_pos) = ztd(idx_is);
                            zhd = this.rec_list(i).work.getAprZhd();
                            this.rec_list(i).work.zwd(idx_pos) = ztd(idx_is) - zhd(idx_pos);
                        end
                        if this.state.flag_tropo_gradient
                            gn = this.ztd_gn(idx_res_av, i);
                            this.rec_list(i).work.tgn(idx_pos) = gn(idx_is);
                            ge = this.ztd_ge(idx_res_av, i);
                            this.rec_list(i).work.tge(idx_pos) = ge(idx_is);
                        end
                        % sigma of the session
                        this.rec_list(i).work.s0 = s0;
                        % residual
                        this.rec_list(i).work.sat.res(:) = 0;
                        this.rec_list(i).work.sat.res(idx_pos, :) = res(idx_is, :, i);
                        
                        this.rec_list(i).work.pushResult();
                        this.rec_list(i).work.updateErrTropo();
                        
                    end
                else
                    this.log.addWarning(sprintf('s0 ( %.4f) too high! not updating the results',s0));
                end
            else
                this.log.addWarning('Not enough receivers (< 2), skipping network solution');
            end
            
        end
    end
end
