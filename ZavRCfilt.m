function dataFlt = ZavRCfilt(data, fc, fsampl, filter_type)
%dataFlt = ZavRCfilt(data, fc, fsampl, filter_type)
%RC-filter. Only highpass
%INPUT
%data - input signal
%fc - cutting frequency (Hz)
%fsampl - sampling frequency (Hz)
%filter_type - type of filter (low - lowpass, high - highpass)
%
%OUTPUT
%dataFlt - output filtered signal
%

%===== righ way (fast) =====%
% Default to allpass if invalid type is selected
b = [1, 0]; a = [1, 0];
z = 1 / tan((2 * pi * fc / fsampl) / 2);

if (strcmp(filter_type, 'high')) %highpass
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The following equations were derived from
    %            s
    % T(s) =  -------
    %          s + 1
    %
    % using Bilinear transform of
    %             1          ( 1 - z^-1 )
    % s -->  -----------  *  ------------
    %         tan(w*T/2)     ( 1 + z^-1 )
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    b(1) = z / (1 + z);
    b(2) = -b(1);
    a(2) = (1 - z) / (1 + z);
else %lowpass %(strcmp(filter_type, 'low'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The following equations were derived from
    %            1
    % T(s) =  -------
    %          s + 1
    %
    % using Bilinear transform of
    %             1          ( 1 - z^-1 )
    % s -->  -----------  *  ------------
    %         tan(w*T/2)     ( 1 + z^-1 )
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    b(1) = 1 / (1 + z);
    b(2) = b(1);
    a(2) = (1 - z) / (1 + z);
end
dataFlt = data;%memory preallocation
for ch = 1:size(data, 2) %run over channels
    for sw = 1:size(data, 3) %run over segments
        dataFlt(:, ch, sw) = filter(b, a, data(:, ch, sw), [], 1);%filtering
    end
end


% %===== short way with same result (slow) =====%
% dataFlt = data;
% RC = 1 / (2 * pi * fc);
% alpha = RC / (RC + (1 / fsampl));
% for ch = 1:size(data, 2)
%     dataFlt(1, ch) = 0;
%     for t = 2:size(data, 1)
%         dataFlt(t, ch) = alpha * (dataFlt(t - 1, ch) + data(t, ch) - data(t - 1, ch));
%     end
% end

