function [filt] = BSfilter(dat,Fs,Fbp,N,type,dir,instabilityfix)

% FT_PREPROC_BANDSTOPFILTER applies a band-stop filter to the data and thereby
% removes the spectral components in the specified frequency band
%
% Use as
%   [filt] = ft_preproc_bandstopfilter(dat, Fsample, Fbp, N, type, dir, instabilityfix)
% where
%   dat        data matrix (Nchans X Ntime)
%   Fsample    sampling frequency in Hz
%   Fbp        frequency band, specified as [Fhp Flp]
%   N          optional filter order, default is 4 (but) or dependent upon
%              frequency band and data length (fir/firls)
%   type       optional filter type, can be
%                'but' Butterworth IIR filter (default)
%                'fir' FIR filter using Matlab fir1 function
%                'firls' FIR filter using Matlab firls function (requires Matlab Signal Processing Toolbox)
%                'brickwall' Frequency-domain filter using Matlab FFT and iFFT function
%   dir        optional filter direction, can be
%                'onepass'         forward filter only
%                'onepass-reverse' reverse filter only, i.e. backward in time
%                'twopass'         zero-phase forward and reverse filter (default)
%                'twopass-reverse' zero-phase reverse and forward filter
%                'twopass-average' average of the twopass and the twopass-reverse
%   instabilityfix optional method to deal with filter instabilities
%                'no'       only detect and give error (default)
%                'reduce'   reduce the filter order
%                'split'    split the filter in two lower-order filters, apply sequentially
%
% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.
%
% Further note that the filter type 'brickwall' filters in the frequency domain,
% but may have severe issues. For instance, it has the implication that the time
% domain signal is periodic. Another issue pertains to that frequencies are
% not well defined over short time intervals; particularly for low frequencies.

% Nyquist frequency
Fn = Fs/2;

% compute filter coefficients
switch type
  case 'but'
    if isempty(N)
      N = 4;
    end
    [B, A] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn], 'stop');
  case 'fir'
    if isempty(N)
      N = 3*fix(Fs / Fbp(1));
      if rem(N,2)==1,   N=N+1;    end
    end
    if N > floor( (size(dat,2) - 1) / 3)
      N=floor(size(dat,2)/3) - 2;
      if rem(N,2)==1,   N=N+1;    end
    end
    [B, A] = fir1(N, [min(Fbp)/Fn max(Fbp)/Fn], 'stop');
  case 'firls' % from NUTMEG's implementation
    if isempty(N)
      N = 3*fix(Fs / Fbp(1));
    end
    if N > floor( (size(dat,2) - 1) / 3)
      N=floor(size(dat,2)/3) - 1;
    end
    f = 0:0.001:1;
    if rem(length(f),2)~=0
      f(end)=[];
    end
    z = ones(1,length(f));
    if(isfinite(min(Fbp)))
      [val,pos1] = min(abs(Fs*f/2 - min(Fbp)));
    else
      [val,pos2] = min(abs(Fs*f/2 - max(Fbp)));
      pos1=pos2;
    end
    if(isfinite(max(Fbp)))
      [val,pos2] = min(abs(Fs*f/2 - max(Fbp)));
    else
      pos2 = length(f);
    end
    z(pos1:pos2) = 0;
    A = 1;
    B = firls(N,f,z); % requires Matlab signal processing toolbox
  case 'brickwall'
    ax = linspace(0, Fs, size(dat,2)); % frequency coefficients
    fl = nearest(ax, min(Fbp))-1; % low cut-off frequency
    fh = nearest(ax, max(Fbp))+1; % high cut-off frequency
    a  = 0; % suppresion rate of frequencies-not-of-interest
    f           = fft(dat,[],2); % FFT
    f(:,fl:fh)  = a.*f(:,fl:fh); % perform band cut-off
    filt        = 2*real(ifft(f,[],2)); % iFFT
    return
  otherwise
    error('unsupported filter type "%s"', type);
end

% demean the data before filtering
meandat = mean(dat,2);
dat = bsxfun(@minus, dat, meandat);

try
  filt = filter_with_correction(B,A,dat,dir);
catch
  switch instabilityfix
    case 'no'
      rethrow(lasterror);
    case 'reduce'
      warning('backtrace', 'off')
      warning('instability detected - reducing the %dth order filter to an %dth order filter', N, N-1);
      warning('backtrace', 'on')
      filt = BSfilter(dat,Fs,Fbp,N-1,type,dir,instabilityfix);
    case 'split'
      N1 = ceil(N/2);
      N2 = floor(N/2);
      warning('backtrace', 'off')
      warning('instability detected - splitting the %dth order filter in a sequential %dth and a %dth order filter', N, N1, N2);
      warning('backtrace', 'on')
      filt1 = BSfilter(dat  ,Fs,Fbp,N1,type,dir,instabilityfix);
      filt  = BSfilter(filt1,Fs,Fbp,N2,type,dir,instabilityfix);
    otherwise
      error('incorrect specification of instabilityfix');
  end % switch
end

% add the mean back to the filtered data
filt = bsxfun(@plus, filt, meandat);


function filt = filter_with_correction(B,A,dat,dir)

% FILTER_WITH_CORRECTION applies a to the data and corrects
% edge-artifacts for one-pass filtering.
%
% Use as
%   [filt] = filter_with_correction(B,A,dat,dir);
% where
%   B,A        filter coefficients
%   dat        data matrix (Nchans X Ntime)
%   dir        optional filter direction, can be
%                'onepass'         forward filter only
%                'onepass-reverse' reverse filter only, i.e. backward in time
%                'twopass'         zero-phase forward and reverse filter (default)
%                'twopass-reverse' zero-phase reverse and forward filter
%                'twopass-average' average of the twopass and the twopass-reverse
%
% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.


poles = roots(A);
if any(abs(poles) >= 1)
  error('Calculated filter coefficients have poles on or outside the unit circle and will not be stable. Try a higher cutoff frequency or a different type/order of filter.');
end

dcGain = sum(B)/sum(A);

[d,N] = size(dat);

switch dir
  case 'onepass'
    offset = dat(:,1);
    dat = dat - repmat(offset,1,N);
    filt = filter(B, A, dat')' + repmat(dcGain*offset, 1, N);
  case 'onepass-reverse'
    offset = dat(:,end);
    dat  = fliplr(dat) - repmat(offset,1,N);
    filt = filter(B, A, dat')';
    filt = fliplr(filt) + repmat(dcGain*offset, 1, N);
  case 'twopass'
    % filtfilt does the correction for us
    filt = filtfilt(B, A, dat')';
  case 'twopass-reverse'
    % filtfilt does the correction for us
    filt = fliplr(filtfilt(B, A, fliplr(dat)')');
  case 'twopass-average'
    % take the average from the twopass and the twopass-reverse
    filt1 = filtfilt(B, A, dat')';
    filt2 = fliplr(filtfilt(B, A, fliplr(dat)')');
    filt  = (filt1 + filt2)/2;
  otherwise
    error('unsupported filter direction "%s"', dir);
end
