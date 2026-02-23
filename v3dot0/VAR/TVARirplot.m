function TVARirplot(IR, TVAR, VARopt, INF, SUP)
%==========================================================================
% Plot regime-specific impulse response functions from a Threshold VAR.
%
% For each shock, one figure is produced. Each subplot shows the IRFs for
% one response variable with BOTH regimes overlaid on the same axes,
% enabling direct cross-regime comparison. If confidence bands are
% provided (INF/SUP from TVARirband), they are drawn as shaded swathes in
% the colours of the respective regimes.
%==========================================================================
% TVARirplot(IR, TVAR, VARopt, INF, SUP)
% -------------------------------------------------------------------------
% INPUT
%   - IR(:,:,:,k) : impulse responses from TVARir
%                   (nsteps x nvar x nvar x nregimes)
%   - TVAR        : TVAR structure (used to read TVARopt.rnames if set in
%                   VARopt, and for the threshold value in the title)
%   - VARopt      : VAR options.  Key fields:
%                     .vnames   - cell array of variable names (required)
%                     .snames   - cell array of shock names
%                     .nsteps   - IRF horizon
%                     .pick     - 0 = all shocks; k = shock k only
%                     .FigSize  - figure size [width, height]
%                     .figname  - file prefix for saving
%                     .quality  - 0=low, 1=high (ghostscript), 2=exportgraphics
%                     .suptitle - 1 = add suptitle
% -------------------------------------------------------------------------
% OPTIONAL INPUT
%   - INF(:,:,:,k) : lower confidence band (from TVARirband), same shape as IR
%   - SUP(:,:,:,k) : upper confidence band (from TVARirband), same shape as IR
% -------------------------------------------------------------------------
% EXAMPLE
%   - See TVARToolbox_Primer.m in "../Primer/"
%==========================================================================
% VAR Toolbox 3.0 - TVAR Extension
%==========================================================================


%% Check inputs
%--------------------------------------------------------------------------
if ~exist('VARopt','var')
    error('TVARirplot: you need to provide VARopt');
end
vnames = VARopt.vnames;
if isempty(vnames)
    error('TVARirplot: add variable names to VARopt.vnames');
end
if isempty(VARopt.snames)
    snames = VARopt.vnames;
else
    snames = VARopt.snames;
end

has_bands = exist('INF','var') && exist('SUP','var') && ~isempty(INF) && ~isempty(SUP);

% Regime labels — read from TVARopt if stored; fall back to generic names
if isfield(VARopt,'rnames') && ~isempty(VARopt.rnames)
    rnames = VARopt.rnames;
else
    rnames = {'Regime 1', 'Regime 2'};
end


%% Retrieve dimensions
%--------------------------------------------------------------------------
[nsteps, nvars, nshocks, nregimes] = size(IR);

pick     = VARopt.pick;
filename = [VARopt.figname 'TVAR_IR_'];
quality  = VARopt.quality;
suptit   = VARopt.suptitle;

if pick < 0 || pick > nvars
    error('TVARirplot: pick value is out of range');
end
shock_start = 1;
shock_end   = nshocks;
if pick > 0
    shock_start = pick;
    shock_end   = pick;
end

% Subplot grid
row = round(sqrt(nvars));
col = ceil(sqrt(nvars));

steps  = 1:nsteps;
x_axis = zeros(1,nsteps);


%% Colours (one per regime)
%--------------------------------------------------------------------------
% cmap(1) = first colour, cmap(2) = second colour, etc.
clr = zeros(nregimes, 3);
for k = 1:nregimes
    clr(k,:) = cmap(k);
end


%% Plot
%--------------------------------------------------------------------------
SwatheOpt        = PlotSwatheOption;
SwatheOpt.trans  = 0.25;   % band transparency
SwatheOpt.marker = '';

for jj = shock_start:shock_end

    FigSize(VARopt.FigSize(1), VARopt.FigSize(2));

    for ii = 1:nvars
        subplot(row, col, ii);

        % --- plot both regimes on same axes ---
        for k = 1:nregimes
            ir_k = IR(:, ii, jj, k);

            % Confidence band
            if has_bands
                inf_k = INF(:, ii, jj, k);
                sup_k = SUP(:, ii, jj, k);
                SwatheOpt.color = clr(k,:);
                PlotSwathe(ir_k, [inf_k sup_k], SwatheOpt); hold on;
            end

            % IRF line
            plot(steps, ir_k, 'Color', clr(k,:), 'LineWidth', 2); hold on;
        end

        % Zero line
        plot(x_axis, '--k', 'LineWidth', 0.5); hold on;

        xlim([1 nsteps]);
        title([vnames{ii} ' to ' snames{jj}], 'FontWeight','bold','FontSize',10);
        set(gca,'Layer','top');

        % Legend on the first subplot only
        if ii == 1
            leg_h = [];
            for k = 1:nregimes
                leg_h(k) = plot(nan, nan, '-', 'Color', clr(k,:), 'LineWidth', 2); %#ok<AGROW>
            end
            legend(leg_h, rnames, 'Location', 'best', 'FontSize', 8);
            legend('boxoff');
        end
    end

    % Optional suptitle
    if suptit
        SupTitle([snames{jj} '  (threshold = ' num2str(TVAR.thresh,'%.3f') ')']);
    end

    % Save figure
    FigName = [filename num2str(jj)];
    if quality == 1
        set(gcf,'Color','w');
        export_fig(FigName, '-pdf', '-painters');
    elseif quality == 2
        exportgraphics(gcf, [FigName '.pdf']);
    elseif quality == 0
        print('-dpdf', '-r100', FigName);
    end

    clf('reset');
end

close all;
