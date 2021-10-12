function dataTrueDC = RecovDCfromPseudo(dataPseudoDC, dFreq)
%
%dataTrueDC = RecovDCfromPseudo(dataPseudoDC, dFreq)
%recovery of DC from pseudo DC signal (Azat)
%
%INPUTS
%dataPseudoDC - pseudoDC signal
%dFreq - descritization frequency (Hz)
%
%OUTPUTS
%dataTrueDC - trueDC signal
%

T = 1 / dFreq;%

Rc = 10;
R = 1;
C = 1;

b = [((Rc + R) * T) + (2 * C * R * Rc), ((Rc + R) * T) - (2 * C * R * Rc)];
a = [(R * T) + (2 * C * R * Rc), (R * T) - (2 * C * R * Rc)];

dataTrueDC = filter(b, a, dataPseudoDC);%trueDC signal recovery
