function [SEGPNTS, DCMP, AVEEGSEGS] = SISegmentation(EEG, Wr, Wd)
% [SEGPNTS, DCMP, AVEEGSEGS] = SISegmentation(EEG, Wr, Wd)
%
% Source-informed segmentation.
% EEG  = (nchan * T * ntrials) EEG matrix with nchan channels, T time samples, and ntrials trials.
% Wr   = reference window size [refwinsize {, slidwinsize, stepsize, overlap}].
% Wd   = decision window size.
% SEGPNTS   = (1 * nsegs) segment boundary indeces (always starting with 1)
% DCMP      = (dcmp.pvals, dcmp.svals, dcmp.svecs, dcmp.eperc) decomposition structure,
%              dcmp.pvals = successful segmentation p-values,
%              dcmp.svals = significant singular values,
%              dcmp.svecs = significant singular vectors, and
%              dcmp.eperc = percentage of energy within the significant space.
% AVEEGSEGS = (nchan * nsegs) average EEG segments.
% If ntrials > 1, outputs are cell arrays of size (ntrials * 1), each cell containing the structures above.
%
% References:
% -----------
% [1] Ali E. Haddad, Laleh Najafizadeh, "Global EEG segmentation using singular value
%     decomposition," 37th Annual International Conference of the IEEE Engineering
%     in Medicine and Biology Society (EMBC), IEEE, 2015, pp. 558-561.
% 
% [2] Ali Haddad, Laleh Najafizadeh, "Multi-scale analysis of the dynamics of brain
%     functional connectivity using EEG," IEEE Biomedical Circuits and Systems
%     Conference (BioCAS), IEEE, 2016, pp. 240-243.
% 
% [3] Ali Haddad, Laleh Najafizadeh, "Source-informed segmentation: Towards capturing
%     the dynamics of brain functional networks through EEG," 50th Asilomar Conference
%     on Signals, Systems and Computers, IEEE, 2016, pp. 1290-1294.
%
% Coded by: Ali Haddad

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
if ntrials > 1
    cellcond  = true;
    SEGPNTS   = cell(ntrials,1);
    DCMP      = cell(ntrials,1);
    AVEEGSEGS = cell(ntrials,1);
else
    cellcond = false;
end
for trial = 1:ntrials
    if dcmpcond
        dcmp.svals = {};
        dcmp.svecs = {};
        dcmp.eperc = {};
        dcmp.pvals = 0;
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
            ep = 1;
            Sr = 0;
            F  = ones(nchan,1) / sqrt(nchan);
        else
            [Ur,Sr] = svd(refwin, 'econ');
            [R,ep]  = curve_elbow(diag(Sr)');
            F       = Ur(:,1:R);
        end
        Er           = sum(((eye(nchan,nchan)-F*F')*refwin).^2);
        Es           = sum(((eye(nchan,nchan)-F*F')*slidwin).^2);
        [dcsn,pvals] = kstest2(Er,Es);
        if dcsn
            dcsnnum     = 1;
            pv          = pvals;
            slidpntminp = slidpnt;
            if dcmpcond
                svals = diag(Sr(1:R,1:R))';
                svecs = F;
                eperc = ep;
            end
            while ((pv<pvals)||(dcsnnum<Wd)||(dcsnnum==1)) && dcsn
                if pv < pvals
                    pvals       = pv;
                    slidpntminp = slidpnt;
                    if dcmpcond
                        svals = diag(Sr(1:R,1:R))';
                        svecs = F;
                        eperc = ep;
                    end
                end
                cursegsize = cursegsize + stepsize;
                slidpnt    = segpnts(segind) + cursegsize - overlap;
                if (segpnts(segind)+cursegsize-1+max(slidwinsize-overlap,0)) <= T
                    refwin  = EEG(:,segpnts(segind):segpnts(segind)+cursegsize-1,trial);
                    slidwin = EEG(:,slidpnt:slidpnt+slidwinsize-1,trial);
                    if sum(abs(refwin(:))) == 0
                        R  = 1;
                        ep = 1;
                        Sr = 0;
                        F  = ones(nchan,1) / sqrt(nchan);
                    else
                        [Ur,Sr] = svd(refwin, 'econ');
                        [R,ep]  = curve_elbow(diag(Sr)');
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
                    dcmp.svals{segind}   = svals;
                    dcmp.svecs{segind}   = svecs;
                    dcmp.eperc{segind}   = eperc;
                    dcmp.pvals(segind+1) = pvals;
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
            dcmp.svals{segind} = 0;
            dcmp.svecs{segind} = ones(nchan,1) / sqrt(nchan);
            dcmp.eperc{segind} = 1;
        else
            [Ur,Sr]            = svd(refwin, 'econ');
            [R,ep]             = curve_elbow(diag(Sr)');
            dcmp.svals{segind} = diag(Sr(1:R,1:R))';
            dcmp.svecs{segind} = Ur(:,1:R);
            dcmp.eperc{segind} = ep;
        end
    end
    if avcond
        augsegpnts = [segpnts,T+1];
        avEEGsegs  = zeros(nchan,segind);
        for seg = 1:segind
            avEEGsegs(:,seg) = mean(EEG(:,augsegpnts(seg):augsegpnts(seg+1)-1,trial),2);
        end
    end
    if cellcond
        SEGPNTS{trial} = segpnts;
        if dcmpcond
            DCMP{trial} = dcmp;
        end
        if avcond
            AVEEGSEGS{trial} = avEEGsegs;
        end
    else
        SEGPNTS = segpnts;
        if dcmpcond
            DCMP = dcmp;
        end
        if avcond
            AVEEGSEGS = avEEGsegs;
        end
    end
end
return