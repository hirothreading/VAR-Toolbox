function TVARirplot(IR, VARopt, INF, SUP)
%==========================================================================
% Plot TVAR impulse responses for two regimes. Two modes:
%   'overlay'  — both regimes overlaid on the same figure (default)
%   'separate' — one figure per regime (follows VARirplot conventions)
%
% IR and bands must be 4-D: (nsteps x nvar x nshocks x nregimes).
% Bands can optionally have a 5th dimension for multiple confidence levels.
%==========================================================================
% TVARirplot(IR, VARopt, INF, SUP)
% -------------------------------------------------------------------------
% INPUT
%   - IR     : (nsteps x nvar x nshocks x nregimes) from TVARir
%   - VARopt : options (vnames, snames, figname, FigSize, etc.)
%              VARopt.tvar_plot_mode: 'overlay' (default) or 'separate'
%              VARopt.regime_names: cell array of regime labels
% -------------------------------------------------------------------------
% OPTIONAL INPUT
%   - INF, SUP : lower/upper bands (same dims as IR, or 5-D with pctg)
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% -------------------------------------------------------------------------


%% Check inputs
%==========================================================================
if ~exist('VARopt','var')
    error('You need to provide VARopt structure');
end
vnames = VARopt.vnames;
if isempty(vnames)
    error('You need to add labels for endogenous variables in VARopt.vnames');
end
if isempty(VARopt.snames)
    snames = VARopt.vnames;
else
    snames = VARopt.snames;
end

% Plot mode
if isfield(VARopt, 'tvar_plot_mode')
    plot_mode = VARopt.tvar_plot_mode;
else
    plot_mode = 'overlay';
end

% Regime names
if isfield(VARopt, 'regime_names') && ~isempty(VARopt.regime_names)
    regime_names = VARopt.regime_names;
else
    regime_names = {'Regime 1', 'Regime 2'};
end


%% Retrieve dimensions
%==========================================================================
[nsteps, nvars, nshocks, nregimes] = size(IR);

filename = [VARopt.figname 'TVAR_IR_'];
quality  = VARopt.quality;
suptitle_flag = VARopt.suptitle;
pick     = VARopt.pick;

if pick < 0 || pick > nvars
    error('The selected shock is non valid');
elseif pick == 0
    pick = 1;
else
    nshocks = pick;
end

% Subplot layout
row = round(sqrt(nvars));
col = ceil(sqrt(nvars));

steps  = 1:nsteps;
x_axis = zeros(1, nsteps);

% Has bands?
has_bands = exist('INF','var') && exist('SUP','var') && ~isempty(INF);


%% Colours
%==========================================================================
% Regime 1: blue
col_r1       = [0, 0.447, 0.741];
swathe_r1_dk = [0.60, 0.75, 0.88];
swathe_r1_lt = [0.80, 0.88, 0.95];

% Regime 2: red
col_r2       = [0.850, 0.325, 0.098];
swathe_r2_dk = [0.92, 0.70, 0.62];
swathe_r2_lt = [0.96, 0.85, 0.80];

SwatheOpt            = PlotSwatheOption;
SwatheOpt.swatheonly = 1;
SwatheOpt.marker     = 'none';


%% Plot: overlay mode
%==========================================================================
if strcmp(plot_mode, 'overlay')
    FigSize(VARopt.FigSize(1), VARopt.FigSize(2));
    for jj = pick:nshocks
        for ii = 1:nvars
            subplot(row, col, ii);

            % --- Regime 1 bands ---
            if has_bands
                plot_bands(IR(:,ii,jj,1), INF, SUP, ii, jj, 1, ...
                    swathe_r1_dk, swathe_r1_lt, SwatheOpt);
            end
            % --- Regime 2 bands ---
            if has_bands
                plot_bands(IR(:,ii,jj,2), INF, SUP, ii, jj, 2, ...
                    swathe_r2_dk, swathe_r2_lt, SwatheOpt);
            end

            % --- Point estimates ---
            plot(steps, IR(:,ii,jj,1), 'LineStyle','-', 'Color',col_r1, 'LineWidth',2); hold on;
            plot(steps, IR(:,ii,jj,2), 'LineStyle','--','Color',col_r2, 'LineWidth',2); hold on;
            plot(x_axis, '--k', 'LineWidth', 0.5); hold on;
            xlim([1 nsteps]);
            title([vnames{ii} ' to ' snames{jj}], 'FontWeight','bold', 'FontSize',10);
            set(gca, 'Layer', 'top');

            if ii == 1
                legend(regime_names{1}, regime_names{2}, 'Location','best');
            end
        end

        % Save
        FigName = [filename num2str(jj)];
        if quality == 1 || quality == 2
            if suptitle_flag == 1
                Alphabet = char('a'+(1:nshocks)-1);
                SupTitle([Alphabet(jj) ') TVAR IR to ' snames{jj}]);
            end
            set(gcf, 'Color', 'w');
            exportgraphics(gcf, [FigName '.pdf'], 'ContentType', 'vector');
        elseif quality == 0
            print('-dpdf', '-r100', FigName);
        end
        clf('reset');
    end
    close all;


%% Plot: separate mode (one figure per regime)
%==========================================================================
elseif strcmp(plot_mode, 'separate')
    for k = 1:nregimes
        if k == 1
            line_col = col_r1;
            sw_dk = swathe_r1_dk;
            sw_lt = swathe_r1_lt;
        else
            line_col = col_r2;
            sw_dk = swathe_r2_dk;
            sw_lt = swathe_r2_lt;
        end

        FigSize(VARopt.FigSize(1), VARopt.FigSize(2));
        for jj = pick:nshocks
            for ii = 1:nvars
                subplot(row, col, ii);

                if has_bands
                    plot_bands(IR(:,ii,jj,k), INF, SUP, ii, jj, k, ...
                        sw_dk, sw_lt, SwatheOpt);
                end
                plot(steps, IR(:,ii,jj,k), 'LineStyle','-', 'Color',line_col, 'LineWidth',2); hold on;
                plot(x_axis, '--k', 'LineWidth', 0.5); hold on;
                xlim([1 nsteps]);
                title([vnames{ii} ' to ' snames{jj}], 'FontWeight','bold', 'FontSize',10);
                set(gca, 'Layer', 'top');
            end

            % Save
            FigName = [filename num2str(jj) '_' num2str(k)];
            if quality == 1 || quality == 2
                if suptitle_flag == 1
                    SupTitle([regime_names{k} ': IR to ' snames{jj}]);
                end
                set(gcf, 'Color', 'w');
                exportgraphics(gcf, [FigName '.pdf'], 'ContentType', 'vector');
            elseif quality == 0
                print('-dpdf', '-r100', FigName);
            end
            clf('reset');
        end
    end
    close all;
end

end


%% ========================================================================
%  HELPER: plot confidence bands for one regime
%  ========================================================================
function plot_bands(IR_line, INF, SUP, ii, jj, k, sw_dk, sw_lt, SwatheOpt)

ndim_inf = ndims(INF);

if ndim_inf == 5
    % Multiple bands: 5th dim is pctg
    nbands = size(INF, 5);
    if nbands == 1
        band_cols = sw_dk;
    else
        band_cols = [linspace(sw_dk(1), sw_lt(1), nbands)', ...
                     linspace(sw_dk(2), sw_lt(2), nbands)', ...
                     linspace(sw_dk(3), sw_lt(3), nbands)'];
    end
    for pp = nbands:-1:1
        SwatheOpt_pp = SwatheOpt;
        SwatheOpt_pp.swathecol = band_cols(pp,:);
        SwatheOpt_pp.linecol   = band_cols(pp,:);
        PlotSwathe(IR_line, [INF(:,ii,jj,k,pp) SUP(:,ii,jj,k,pp)], SwatheOpt_pp);
        hold on;
    end
else
    % Single band
    SwatheOpt_s = SwatheOpt;
    SwatheOpt_s.swathecol = sw_dk;
    SwatheOpt_s.linecol   = sw_dk;
    PlotSwathe(IR_line, [INF(:,ii,jj,k) SUP(:,ii,jj,k)], SwatheOpt_s);
    hold on;
end

end
