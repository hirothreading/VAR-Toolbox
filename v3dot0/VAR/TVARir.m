function [IR, TVAR] = TVARir(TVAR, VARopt)
%==========================================================================
% Compute regime-specific impulse response functions (IRFs) for a
% Threshold VAR estimated with TVARmodel. Four identification schemes are
% supported: zero contemporaneous restrictions ('short'), zero long-run
% restrictions ('long'), sign restrictions ('sign'), and external
% instruments / proxy SVAR ('iv').
%
% For 'short', 'long', and 'sign', this function delegates to the existing
% VARir for each regime. For 'iv' (Proxy TVAR), the instrument is first
% trimmed to regime-specific observations before the two-stage estimation.
%
% The IRFs for the two regimes are stacked in a 4-D output array so that
% they can be passed directly to TVARirplot and TVARirband.
%==========================================================================
% [IR, TVAR] = TVARir(TVAR, VARopt)
% -------------------------------------------------------------------------
% INPUT
%   - TVAR   : structure, result of TVARmodel.m
%   - VARopt : VAR options (see VARoption / TVARoption).
%              VARopt.ident controls identification:
%                'short' - zero contemporaneous restrictions (Cholesky)
%                'long'  - zero long-run restrictions
%                'sign'  - sign restrictions (TVAR.regime{k}.B must be
%                          pre-populated via SR before calling TVARir)
%                'iv'    - external instrument; requires TVAR.IV to be set
% -------------------------------------------------------------------------
% OUTPUT
%   - IR(:,:,:,k) : impulse responses for regime k
%                   Dimensions: (nsteps x nvar x nvar x nregimes)
%                   IR(:,:,:,1) = regime 1,  IR(:,:,:,2) = regime 2
%   - TVAR        : updated structure; TVAR.regime{k}.B, .Biv, .PSI, .Fp
%                   are populated for each regime
% -------------------------------------------------------------------------
% NOTE ON 'sign' IDENTIFICATION
%   Sign restrictions require TVAR.regime{k}.B to be pre-computed for each
%   regime. Call SR(TVAR.regime{k}, SIGN, VARopt) for k=1,2 and store the
%   result in TVAR.regime{k}.B before calling TVARir.
% -------------------------------------------------------------------------
% EXAMPLE
%   - See TVARToolbox_Primer.m in "../Primer/"
%==========================================================================
% VAR Toolbox 3.0 - TVAR Extension
%==========================================================================


%% Check inputs
%--------------------------------------------------------------------------
if ~exist('TVAR','var')
    error('You need to provide TVAR structure, result of TVARmodel');
end
if ~exist('VARopt','var')
    error('You need to provide VARopt structure');
end

ident    = VARopt.ident;
nvar     = TVAR.nvar;
nsteps   = VARopt.nsteps;
nregimes = TVAR.nregimes;

% Check instrument is available for iv identification
if strcmp(ident, 'iv')
    if isempty(TVAR.IV)
        error('TVARir: VARopt.ident=''iv'' requires TVAR.IV to be set');
    end
    % TVAR.IV should have nobs rows (aligned with the original ENDO before
    % lag-trimming), matching the convention in VARir.
    if size(TVAR.IV, 1) ~= size(TVAR.ENDO, 1)
        error('TVARir: TVAR.IV must have the same number of rows as TVAR.ENDO (%d)', size(TVAR.ENDO,1));
    end
end

% Pre-allocate output
IR = nan(nsteps, nvar, nvar, nregimes);


%% Loop over regimes
%--------------------------------------------------------------------------
for k = 1:nregimes

    reg = TVAR.regime{k};

    % ------------------------------------------------------------------
    % Identification with external instrument ('iv') — handled here to
    % avoid the nlag-offset ambiguity present in VARir
    % ------------------------------------------------------------------
    if strcmp(ident, 'iv')

        % Regime-specific instrument: rows that fall in this regime.
        % TVAR.regime{k}.obs_idx contains the row indices into Y/X (which
        % have already had the first nlag rows of ENDO removed).
        % The full instrument TVAR.IV has nobs rows (same as ENDO).
        % To align: row t of Y corresponds to row nlag+t of ENDO, so the
        % instrument value for that observation is IV(nlag+t).
        % Therefore IV rows for this regime = IV(nlag + obs_idx).
        nlag  = TVAR.nlag;
        IV_k  = TVAR.IV(nlag + reg.obs_idx, :);

        up = reg.resid(:, 1);       % residuals to be instrumented (nobs_k x 1)
        uq = reg.resid(:, 2:end);   % other residuals               (nobs_k x (nvar-1))

        % Align IV with residuals (CommonSample handles leading/trailing NaNs)
        [aux, ~, ~] = CommonSample([up IV_k]);
        p  = aux(:, 1);
        q  = uq(end-length(p)+1:end, :);
        pq = [p q];
        Z  = aux(:, 2:end);

        if isempty(Z) || size(Z,1) < 2
            error('TVARir: insufficient non-NaN instrument observations in regime %d for iv identification', k);
        end

        % First stage: regress instrumented residual on instrument
        FirstStage = OLSmodel(p, Z);
        p_hat = FirstStage.yhat;

        % Second stage: recover remaining columns of B
        Biv_k    = zeros(nvar, 1);
        Biv_k(1) = 1;   % normalisation
        sqsp     = zeros(nvar-1, 1);
        for ii = 2:nvar
            SecondStage = OLSmodel(q(:, ii-1), p_hat);
            Biv_k(ii)   = SecondStage.beta(2);
            sqsp(ii-1)  = SecondStage.beta(2);
        end

        % Scale shock size (Gertler & Karadi 2015 fn 4 / Mertens & Ravn 2013)
        sigma_b = (1/(length(pq) - reg.ntotcoeff)) * ...
            (pq - repmat(mean(pq), size(pq,1), 1))' * ...
            (pq - repmat(mean(pq), size(pq,1), 1));
        s21s11  = sqsp;
        S11     = sigma_b(1, 1);
        S21     = sigma_b(2:end, 1);
        S22     = sigma_b(2:end, 2:end);
        Q       = s21s11*S11*s21s11' - (S21*s21s11' + s21s11*S21') + S22;
        sp      = sqrt(S11 - (S21 - s21s11*S11)' * (Q \ (S21 - s21s11*S11)));
        Biv_k   = Biv_k * sp;

        % Store in regime struct and set B (first column only; rest are zeros)
        reg.Biv    = Biv_k;
        reg.B      = zeros(nvar, nvar);
        reg.B(:,1) = Biv_k;
        reg.IV     = IV_k;
        reg.FirstStage = FirstStage;
        reg.sigma_b    = sigma_b;

        % Delegate to VARir using 'short' stub — B is pre-set, so we
        % temporarily switch ident to 'sign' (which trusts VAR.B directly)
        VARopt_k        = VARopt;
        VARopt_k.ident  = 'sign';   % VARir 'sign' path: uses pre-set VAR.B
        [IR_k, reg]     = VARir(reg, VARopt_k);

    % ------------------------------------------------------------------
    % All other identification schemes — delegate directly to VARir
    % ------------------------------------------------------------------
    else
        reg.IV = [];
        [IR_k, reg] = VARir(reg, VARopt);
    end

    IR(:,:,:,k)  = IR_k;
    TVAR.regime{k} = reg;
end
