function report = runP0Audit(mode)
% RUNP0AUDIT Top-level entry point for P0 8-point audit.
% Usage:
%   runP0Audit          % runs fast suite (default)
%   runP0Audit('fast')  % runs fast suite (N=20, M=100, NSteps=90, t0=31, D=20/40)
%   runP0Audit('full')  % runs full single-seed audit (N=20, M=100, NSteps=2200, t0=1401, D=20/120)

if nargin < 1 || isempty(mode)
    mode = 'fast';
end

switch lower(mode)
    case 'fast'
        report = runP0AuditFast();
    case 'full'
        report = runP0AuditFullSingleSeed();
    otherwise
        error('Unknown audit mode: %s. Use ''fast'' or ''full''.', mode);
end
end
