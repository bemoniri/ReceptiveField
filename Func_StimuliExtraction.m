function [out]  = Func_StimuliExtraction( E, stimuli, freq )
% StimuliExtraction_Func
%  Input : events = spike times     and      stimuli = msq1D.mat file
%  Output : A 16*16*N matrix of stimuli that caused a spike ( N = total
%  number of spikes)
events= ceil(E .*0.0001.*freq);

j = 0;
for i = 1 : length(E)
    if(events(i) > 16)
        j = j + 1;
        out(:, :, j)  = stimuli(events(i) - 15 : events(i), :);
    end
end
end
