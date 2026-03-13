function SaveFigure(path,~,type)
% =======================================================================
% Saves figure to specified path using exportgraphics
% =======================================================================
% SaveFigure(path,~,type)
% -----------------------------------------------------------------------
% INPUT
%   - path: path where to save the file, without extension [char]
% -----------------------------------------------------------------------
% OPTIONAL INPUT
%   - type: pdf, png, eps [dflt=pdf] [char]
% =======================================================================
% VAR Toolbox 3.0
% Ambrogio Cesa-Bianchi
% ambrogiocesabianchi@gmail.com
% March 2015. Updated March 2026
% -----------------------------------------------------------------------

if ~exist('type','var')
    type = 'pdf';
end

set(gcf, 'Color', 'w');
exportgraphics(gcf, [path '.' type], 'ContentType', 'vector');
