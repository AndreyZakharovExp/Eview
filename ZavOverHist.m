function [count, centbins] = ZavOverHist(data, binStep, segmEdge)
%[count, centbins] = ZavOverHist(data, binStep, segmEdge)
%calculating histogram with overlapped bins
%
%INPUTS
%data - sourse data to be treated (time moments, ms)
%binStep - binning parameters; binStp(1) - bin size (ms), bisStp(2) - precision factor (must be integer) 
%segmEdge - time range (ms or same with data). Begin- and end-time to be calculated within
%
%OUTPUTS
%count - calculated count for data
%centbins - resultant binning (centres, ms as in data)
%

binStep(2) = round(binStep(2));%rounding (precision factor must be integer)
bnstp2 = (binStep(1) / 2);%half of base bin
prcbin = binStep(1) / binStep(2);%precise bin

basbins = segmEdge(1):binStep(1):(segmEdge(2) + binStep(1));%base binning - edges
centbins = segmEdge(1):prcbin:segmEdge(2);%base binning - centers
%[centbins(1) - segmEdge(1), centbins(end) - segmEdge(2)] => [left lag, right lag]

% %= old verstion of centbins =%
% centbins = (basbins(1) + bnstp2):prcbin:(basbins(end) - bnstp2 + (prcbin * (binStep(2) - 1)));%base binning - centers
% centbins = centbins - mean([centbins(1) - segmEdge(1), centbins(end) - segmEdge(2)]);
% %= end of old verstion of centbins =%

basbins = basbins - basbins(1);%shifted base binning
cntLen = length(centbins);%number of binning centers
count = zeros(cntLen, 1);%density data
for k = 1:binStep(2) %run over precision steps
    bfr = histc(data, (centbins(k) - bnstp2) + basbins);%source histogram
    jj = k:binStep(2):cntLen;%indexes for current histc
    count(jj) = bfr(1:length(jj));%excepting last bins because of histc count any values that match right edges
end

