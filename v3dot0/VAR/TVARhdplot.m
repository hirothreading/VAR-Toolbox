function TVARhdplot(HD, VARopt, TVAR)
%==========================================================================
% Plot the historical decomposition from TVARhd. Optionally shades
% background to indicate regime membership.
%==========================================================================
% TVARhdplot(HD, VARopt, TVAR)
% -------------------------------------------------------------------------
% INPUT
%   - HD     : structure from TVARhd
%   - VARopt : options (vnames, snames, figname, FigSize, etc.)
%   - TVAR   : (optional) TVAR structure; if provided, regime shading is
%              added to each subplot
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% -------------------------------------------------------------------------


%% Check inputs
%==========================================================================
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
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


%% Retrieve parameters
%==========================================================================
filename = [VARopt.figname 'TVAR_HD_'];
quality  = VARopt.quality;
suptitle = VARopt.suptitle;
pick     = VARopt.pick;

[nsteps, nvars, nshocks] = size(HD.shock);

% If one variable is chosen
if pick < 0 || pick > nvars
    error('The selected variable is non valid');
else
    if pick == 0
        pick = 1;
    else
        nvars = pick;
    end
end

% Regime shading
do_shade = exist('TVAR','var') && ~isempty(TVAR);


%% Plot
%==========================================================================
for ii = pick:nvars
    FigSize(VARopt.FigSize(1), VARopt.FigSize(2));

    % Regime shading (grey bands for regime 2)
    if do_shade
        add_regime_shading(TVAR, nsteps);
    end

    H = AreaPlot(squeeze(HD.shock(:,:,ii))); hold on;
    h = plot(sum(squeeze(HD.shock(:,:,ii)),2), '-k', 'LineWidth', 2);
    if ~isempty(VARopt.firstdate)
        DatesPlot(VARopt.firstdate, nsteps, 8, VARopt.frequency);
    end
    xlim([1 nsteps]);
    set(gca, 'Layer', 'top');
    title(vnames{ii}, 'FontWeight', 'bold', 'FontSize', 10);

    % Save
    FigName = [filename num2str(ii)];
    if quality
        if suptitle == 1
            Alphabet = char('a'+(1:nvars)-1);
            SupTitle([Alphabet(ii) ') TVAR HD of ' vnames{ii}]);
        end
        opt = LegOption; opt.handle = [H(1,:) h];
        LegSubplot([snames(:)' {'Data'}], opt);
        set(gcf, 'Color', 'w');
        exportgraphics(gcf, [FigName '.pdf'], 'ContentType', 'vector');
    else
        legend([H(1,:) h], [snames(:)' {'Data'}]);
        print('-dpdf', '-r100', FigName);
    end
    clf('reset');
end

close all;
end


%% ========================================================================
%  HELPER: add grey shading for regime 2 periods
%  ========================================================================
function add_regime_shading(TVAR, nsteps)
    nlag = TVAR.nlag;
    idx1 = TVAR.regime_idx;   % logical (nobse x 1), true = regime 1
    % HD has nsteps = nobse + nlag rows. First nlag rows are NaN.
    % idx1 corresponds to rows (nlag+1 : nsteps) of HD.
    % Build full index with NaN-padded head
    full_idx = true(nsteps, 1);
    nobse = length(idx1);
    full_idx(nlag+1:nlag+nobse) = idx1;

    % Shade regime 2 intervals (where full_idx == false)
    hold on;
    ylims = get(gca, 'YLim');
    t = 1;
    while t <= nsteps
        if ~full_idx(t)
            t_start = t;
            while t <= nsteps && ~full_idx(t)
                t = t + 1;
            end
            patch([t_start-0.5 t-0.5 t-0.5 t_start-0.5], ...
                  [ylims(1) ylims(1) ylims(2) ylims(2)], ...
                  [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        else
            t = t + 1;
        end
    end
    hold on;
end
