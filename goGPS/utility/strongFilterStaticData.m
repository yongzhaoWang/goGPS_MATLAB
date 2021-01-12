function [data, lid_ko, trend, spline] = strongFilterStaticData(data, robustness_perc, n_sigma, spline_base)
% Returns the data removing outliers (spikes)
%
% INPUT:
%   data                column array of values
%   robustness_perc     maximum percentage of date with no outliers
%
% SYNTAX:
%   [data, id_ko] = strongFilterStaticData(data, robustness_perc)
    
%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b8
%
%--------------------------------------------------------------------------
%  Copyright (C) 2020 Gatti Andrea, Giulio Tagliaferro, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti ...
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

    if any(data(:,end))
        if nargin < 2
            robustness_perc = 0.8;
        end
        if nargin < 3
            n_sigma = 6;
        end
        if nargin < 4 || numel(spline_base) ~= 2
            spline_base = [28, 7, 3.5];
        end
        flag_time = false;
        idf = [];
        if size(data,2) >= 2
            if size(data,2) >= 3
                data_var = data(:,3);
            end
            % Suppose regularly sampled data, fill missing epochs with nan
            time = data(:,1);
            rate = round(median(diff(time*86400)))/86400;
            time_full = linspace(time(1), time(end), round((time(end) - time(1)) / rate + 1))';
            [~, idf, idr] = intersect(round((time_full-rate/2)/rate), round((time-rate/2)/rate));
            tmp = data(:,2);
            data = nan(numel(time_full), 1);
            if numel(idr) < numel(tmp)
                Core.getLogger.addWarning('StrongFilter is loosing some observations out of sync');
            end
            data(idf) = tmp(idr);
            flag_time = 1;
        end

        % Compute a trend "robust" using the robustness_perc of data
        [tmp, trend] = strongDeTrend(data, robustness_perc, 1-((1-robustness_perc)/2), n_sigma);

        if any(tmp) && flag_time && (numel(data(idf)) > 4)
            if (numel(tmp(idf)) > 11)
                spline_base = min(floor(time(end)-time(1)), spline_base);
                warning off;
                % Perform a bit of outlier detection before computing splines
                thr = 6 * perc(abs(tmp), 0.8);
                lid_ok = abs(tmp(idf)) < thr;

                % Computer long splines (reduce the signal, montly splines)
                if size(data,2) >= 3
                    [~, ~, ~, long_spline] = splinerMat(time(lid_ok), [data(idf(lid_ok)) data_var(lid_ok)], spline_base(1), 1e-5, time); % long splines
                else
                    [~, ~, ~, long_spline] = splinerMat(time(lid_ok), [data(idf(lid_ok)) tmp(idf(lid_ok)).^2], spline_base(1), 1e-5, time); % long splines
                end

                % Keep in tmp the reduced value
                tmp(idf) = data(idf) - long_spline;

                % Compute medium splines (reduce the signal weekly splines)
                [~, ~, ~, spline] = splinerMat(time(lid_ok), [tmp(idf(lid_ok)) abs(tmp(idf(lid_ok)))], spline_base(2), 1e-5, time); % medium splines
                [~, ~, ~, spline] = splinerMat(time(lid_ok), [tmp(idf(lid_ok)) abs(tmp(idf(lid_ok)) - spline(lid_ok))], spline_base(2), 1e-5, time); % medium splines

                % These are the medium long frequencies, I reduce the signal so that the interpolation will be more stable
                long_spline = long_spline + spline;

                % Keep in tmp the reduced value
                tmp(idf) = data(idf) - long_spline;

                % Remove high frequencies
                thr = 12 * strongStd(tmp(idf), robustness_perc);
                lid_ok = abs(tmp(idf)) < thr;
                if sum(lid_ok) > 2
                    [~, ~, ~, spline] = splinerMat(time(lid_ok), [tmp(idf(lid_ok)) data_var(lid_ok)], spline_base(3), 1e-2, time); % short splines
                end
                warning on;

                spline = spline + long_spline;
                tmp = data(idf) - spline;
            else
                tmp = tmp(idf);
                spline = trend(idf);
            end
        else
            if flag_time
                spline = trend(idf);
            else
                spline = trend;
            end
        end

        if flag_time
            data = data(idf);
            trend = trend(idf);
        end

        % Outlier detection based on the interpolation
        thr = n_sigma * strongStd(tmp, robustness_perc);
        lid_ko = abs(tmp) > thr;

        % figure; plot(data, 'Color', [0.5 0.5 0.5]);
        data(lid_ko) = nan;
        % hold on; plot(data, '.-b', 'LineWidth', 2)
        % plot(tmp,'g');
    else
        % no data
        data = data(:,end);
        lid_ko = true(size(data));
        trend = data;
        spline = data;
    end
end