function [SEGPNTS, DCMP, AVEEGSEGS] = SISegmentation(EEG, Wr, Wd)
% [SEGPNTS, DCMP, AVEEGSEGS] = SISegmentation(EEG, Wr, Wd)
%
% Source-informed segmentation.
% EEG  = (nchan * T * ntrials) EEG matrix with nchan channels, T time samples, and ntrials trials.
% Wr   = reference window size [refwinsize {, slidwinsize, stepsize, overlap}].
% Wd   = (optional, default Wd = 1) decision window size.
% SEGPNTS   = (ntrials * 1) cell array with (1 * nsegs) segment boundary indeces (always starting with 1).
% DCMP      = (DCMP.svals, DCMP.svecs, DCMP.erat, DCMP.pval) decomposition structure,
%              DCMP.svals = (ntrials * 1) cell array of (1 * nsegs) cell arrays with (1 * R) significant singular values,
%              DCMP.svecs = (ntrials * 1) cell array of (1 * nsegs) cell arrays with (nchan * R) significant singular vectors,
%              DCMP.erat  = (ntrials * 1) cell array with (1 * nsegs) ratios of energy within the significant space,
%              DCMP.pval  = (ntrials * 1) cell array with (1 * nsegs) successful segmentation p-values.
% AVEEGSEGS = (ntrials * 1) cell array of (1 * nsegs) cell arrays of (nchan * 1) average EEG segments.

[nchan,T,ntrials] = size(EEG);
refwinsize        = Wr(1);
winopt            = length(Wr);
if winopt < 2
    slidwinsize = refwinsize;
else
    slidwinsize = Wr(2);
end
if winopt < 3
    stepsize = 1;
else
    stepsize = Wr(3);
end
if winopt < 4
    overlap = 0;
else
    overlap = Wr(4);
end
if (overlap>=slidwinsize) || (overlap>=refwinsize)
    error('overlap must be less than refwinsize and slidwinsize!');
end
if nargin < 3
    Wd = 1;
end
if nargout > 1
    dcmpcond = true;
else
    dcmpcond = false;
end
if nargout > 2
    avcond = true;
else
    avcond = false;
end
SEGPNTS = cell(ntrials,1);
if dcmpcond
    DCMP.svals = cell(ntrials,1);
    DCMP.svecs = cell(ntrials,1);
    DCMP.erat  = cell(ntrials,1);
    DCMP.pval  = cell(ntrials,1);
end
if avcond
    AVEEGSEGS = cell(ntrials,1);
end
for trial = 1:ntrials
    if dcmpcond
        DCMP.svals{trial} = {};
        DCMP.svecs{trial} = {};
        DCMP.erat{trial}  = [];
        DCMP.pval{trial}  = 0;
    end
    segind     = 1;
    segpnts    = 1;
    cursegsize = refwinsize;
    while (segpnts(segind)+cursegsize-1+max(slidwinsize-overlap,0)) <= T
        refwin  = EEG(:,segpnts(segind):segpnts(segind)+cursegsize-1,trial);
        slidpnt = segpnts(segind) + cursegsize - overlap;
        slidwin = EEG(:,slidpnt:slidpnt+slidwinsize-1,trial);
        if sum(abs(refwin(:))) == 0
            R  = 1;
            er = 1;
            Sr = 0;
            F  = ones(nchan,1) / sqrt(nchan);
        else
            [Ur,Sr] = svd(refwin, 'econ');
            [R,er]  = curve_elbow(diag(Sr)');
            F       = Ur(:,1:R);
        end
        Er          = sum(((eye(nchan,nchan)-F*F')*refwin).^2);
        Es          = sum(((eye(nchan,nchan)-F*F')*slidwin).^2);
        [dcsn,pval] = kstest2(Er,Es);
        if dcsn
            dcsnnum     = 1;
            pv          = pval;
            slidpntminp = slidpnt;
            if dcmpcond
                svals = diag(Sr(1:R,1:R))';
                svecs = F;
                erat  = er;
            end
            while ((pv<pval)||(dcsnnum<Wd)||(dcsnnum==1)) && dcsn
                if pv < pval
                    pval        = pv;
                    slidpntminp = slidpnt;
                    if dcmpcond
                        svals = diag(Sr(1:R,1:R))';
                        svecs = F;
                        erat  = er;
                    end
                end
                cursegsize = cursegsize + stepsize;
                slidpnt    = segpnts(segind) + cursegsize - overlap;
                if (segpnts(segind)+cursegsize-1+max(slidwinsize-overlap,0)) <= T
                    refwin  = EEG(:,segpnts(segind):segpnts(segind)+cursegsize-1,trial);
                    slidwin = EEG(:,slidpnt:slidpnt+slidwinsize-1,trial);
                    if sum(abs(refwin(:))) == 0
                        R  = 1;
                        er = 1;
                        Sr = 0;
                        F  = ones(nchan,1) / sqrt(nchan);
                    else
                        [Ur,Sr] = svd(refwin, 'econ');
                        [R,er]  = curve_elbow(diag(Sr)');
                        F       = Ur(:,1:R);
                    end
                    Er        = sum(((eye(nchan,nchan)-F*F')*refwin).^2);
                    Es        = sum(((eye(nchan,nchan)-F*F')*slidwin).^2);
                    [dcsn,pv] = kstest2(Er,Es);
                    if dcsn
                        dcsnnum = dcsnnum + 1;
                    end
                else
                    break
                end
            end
            if dcsnnum >= Wd
                if dcmpcond
                    DCMP.svals{trial}{segind}  = svals;
                    DCMP.svecs{trial}{segind}  = svecs;
                    DCMP.erat{trial}(segind)   = erat;
                    DCMP.pval{trial}(segind+1) = pval;
                end
                segind          = segind + 1;
                segpnts(segind) = slidpntminp; %#ok<AGROW>
                cursegsize      = refwinsize;
            else
                cursegsize = cursegsize + stepsize;
            end
        else
            cursegsize = cursegsize + stepsize;
        end
    end
    if dcmpcond
        refwin = EEG(:,segpnts(segind):T,trial);
        if sum(abs(refwin(:))) == 0
            DCMP.svals{trial}{segind} = 0;
            DCMP.svecs{trial}{segind} = ones(nchan,1) / sqrt(nchan);
            DCMP.erat{trial}(segind)  = 1;
        else
            [Ur,Sr]                   = svd(refwin, 'econ');
            [R,er]                    = curve_elbow(diag(Sr)');
            DCMP.svals{trial}{segind} = diag(Sr(1:R,1:R))';
            DCMP.svecs{trial}{segind} = Ur(:,1:R);
            DCMP.erat{trial}(segind)  = er;
        end
    end
    if avcond
        augsegpnts       = [segpnts,T+1];
        AVEEGSEGS{trial} = cell(1,segind);
        for seg = 1:segind
            AVEEGSEGS{trial}{seg} = mean(EEG(:,augsegpnts(seg):augsegpnts(seg+1)-1,trial),2);
        end
    end
    SEGPNTS{trial} = segpnts;
end
return