function TVARvdplot(VD, VARopt)
%==========================================================================
% Plot TVAR variance decompositions per regime using stacked area plots.
%==========================================================================
% TVARvdplot(VD, VARopt)
% -------------------------------------------------------------------------
% INPUT
%   - VD     : (nsteps x nshocks x nvar x nregimes) from TVARvd
%   - VARopt : options (vnames, snames, figname, FigSize, etc.)
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% -------------------------------------------------------------------------

vnames = VARopt.vnames;
if isempty(vnames)
    error('You need to add labels for endogenous variables in VARopt.vnames');
end
if isempty(VARopt.snames)
    snames = vnames;
else
    snames = VARopt.snames;
end

if isfield(VARopt, 'regime_names') && ~isempty(VARopt.regime_names)
    regime_names = VARopt.regime_names;
else
    regime_names = {'Regime 1', 'Regime 2'};
end

[nsteps, nshocks, nvars, nregimes] = size(VD);
filename = [VARopt.figname 'TVAR_VD_'];
quality  = VARopt.quality;

% Subplot layout
row = round(sqrt(nvars));
col = ceil(sqrt(nvars));

for k = 1:nregimes
    FigSize(VARopt.FigSize(1), VARopt.FigSize(2));
    for ii = 1:nvars
        subplot(row, col, ii);
        H = AreaPlot(VD(:,:,ii,k));
        xlim([1 nsteps]); ylim([0 100]);
        title([vnames{ii} ' — ' regime_names{k}], 'FontWeight','bold', 'FontSize',10);
        set(gca, 'Layer', 'top');
    end
    FigName = [filename num2str(k)];
    if quality
        opt = LegOption; opt.handle = H(1,:);
        LegSubplot(snames, opt);
        set(gcf, 'Color', 'w');
        exportgraphics(gcf, [FigName '.pdf'], 'ContentType', 'vector');
    else
        legend(H(1,:), snames);
        print('-dpdf', '-r100', FigName);
    end
    clf('reset');
end
close all;
