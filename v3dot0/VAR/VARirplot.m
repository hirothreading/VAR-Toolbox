function VARirplot(IR,VARopt,INF,SUP)
% =======================================================================
% Plot the IRs computed with VARir
% =======================================================================
% VARirplot(IR,VARopt,vnames,INF,SUP)
% -----------------------------------------------------------------------
% INPUT
%   - IR(:,:,:) : matrix with IRF (H horizons, N variables, N shocks)
%   - VARopt: options of the VAR (see VARopt from VARmodel)
% -----------------------------------------------------------------------
% OPTIONAL INPUT
%   - INF: lower error band
%   - SUP: upper error band
% -----------------------------------------------------------------------
% EXAMPLE
%   - See VARToolbox_Code.m in "../Primer/"
% =======================================================================
% VAR Toolbox 3.0
% Ambrogio Cesa-Bianchi
% ambrogiocesabianchi@gmail.com
% March 2012. Updated April 2021
% -----------------------------------------------------------------------


%% Check inputs
%================================================
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end
% If there is VARopt get the vnames
vnames = VARopt.vnames;
% Check they are not empty
if isempty(vnames)
    error('You need to add label for endogenous variables in VARopt');
end
% Define shock names
if isempty(VARopt.snames)
    snames = VARopt.vnames;
else
    snames = VARopt.snames;
end

%% Retrieve and initialize variables 
%================================================
filename = [VARopt.figname 'IR_'];
quality = VARopt.quality;
suptitle = VARopt.suptitle;
pick = VARopt.pick;

% Initialize IR matrix
[nsteps, nvars, nshocks] = size(IR);

% If one shock is chosen, set the right value for nshocks
if pick<0 || pick>nvars
    error('The selected shock is non valid')
else
    if pick==0
        pick=1;
    else
        nshocks = pick;
    end
end

% Define the rows and columns for the subplots
row = round(sqrt(nvars));
col = ceil(sqrt(nvars));

% Define a timeline
steps = 1:1:nsteps;
x_axis = zeros(1,nsteps);


%% Plot
%================================================
% Colours
blue       = [0, 0.447, 0.741];   % MATLAB default blue for point estimate
grey_dark  = [0.50, 0.50, 0.50];  % grey for narrowest band
grey_light = [0.80, 0.80, 0.80];  % light grey for widest band

SwatheOpt            = PlotSwatheOption;
SwatheOpt.swatheonly = 1;   % point estimate plotted separately below
SwatheOpt.marker     = 'none';
FigSize(VARopt.FigSize(1),VARopt.FigSize(2))
for jj=pick:nshocks
    for ii=1:nvars
        subplot(row,col,ii);
        if exist('INF','var') && exist('SUP','var')
            if ndims(INF) == 4
                % Multiple confidence bands (vector pctg):
                % pctg stored ascending, so band pp=nbands is widest (lightest grey)
                nbands = size(INF, 4);
                if nbands == 1
                    band_cols = grey_dark;
                else
                    band_cols = [linspace(grey_dark(1), grey_light(1), nbands)', ...
                                 linspace(grey_dark(2), grey_light(2), nbands)', ...
                                 linspace(grey_dark(3), grey_light(3), nbands)'];
                end
                % Plot widest band first (background), narrowest last (foreground)
                for pp = nbands:-1:1
                    SwatheOpt_pp           = SwatheOpt;
                    SwatheOpt_pp.swathecol = band_cols(pp,:);
                    SwatheOpt_pp.linecol   = band_cols(pp,:);
                    PlotSwathe(IR(:,ii,jj), [INF(:,ii,jj,pp) SUP(:,ii,jj,pp)], SwatheOpt_pp); hold on;
                end
            else
                % Single band (scalar pctg)
                SwatheOpt_s           = SwatheOpt;
                SwatheOpt_s.swathecol = grey_dark;
                SwatheOpt_s.linecol   = grey_dark;
                PlotSwathe(IR(:,ii,jj), [INF(:,ii,jj) SUP(:,ii,jj)], SwatheOpt_s); hold on;
            end
        end
        plot(steps, IR(:,ii,jj), 'LineStyle','-','Color',blue,'LineWidth',2); hold on
        plot(x_axis,'--k','LineWidth',0.5); hold on
        xlim([1 nsteps]);
        title([vnames{ii} ' to ' snames{jj}], 'FontWeight','bold','FontSize',10);
        set(gca, 'Layer', 'top');
    end
    % Save
    FigName = [filename num2str(jj)];
    if quality==1
        if suptitle==1
            Alphabet = char('a'+(1:nshocks)-1);
            SupTitle([Alphabet(jj) ') IR to a shock to '  vnames{jj}])
        end
        set(gcf, 'Color', 'w');
        export_fig(FigName,'-pdf','-painters')
    elseif quality==2
        exportgraphics(gcf,[FigName '.pdf'])
    elseif quality==0
        print('-dpdf','-r100',FigName);
    end
    clf('reset');
end

close all