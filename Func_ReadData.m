function [output, freq]  = Func_ReadData( neuronCode )
%FUNC_READDATA
%  Input : Neuron code as a string called neuronCode
%  Output : A struct containing spike times and sa0 headers related to
%  msq1D stimulus

% output(i).events : Spike times related to i th experiment of the neuron
% output(i).header : Header of i th experiment's sa0 file of the  neuron
% freq             : Stimuli Frequency

output = struct('events',{},'hdr',{},'rate', {});
frames = 32767;

% Searching for the neuron
cd Data/Spike_and_Log_Files
arg = ['*', neuronCode, '*'];
folder = dir (arg);
if(length(folder)>1)
    cd ../.. %Go up two directories
   error('More than one neuron found with this criteria'); 
end

cd (folder.name);
addpath('../../../MatlabFunctions/tview','../../../MatlabFunctions/fileload', '../../../Data/Stimulus_Files','../../../Data/Spike_and_Log_Files', '.');

try
    content = dir;
    j = 1;
    for i = 1 : length(content)
       if(~isempty(strfind(content(i).name, 'msq1D.log')) || ~isempty(strfind(content(i).name, 'msq1d.log')))
            ce = importdata(content(i).name);
            u = ce.textdata;
            freq_str = u{13}(31:end);
            f(j) = str2double(freq_str);
       end
    end

    freq = mean(f);

    total_time = frames/freq;
    exp_index = 1;

    % Opening the related msq1d.sa0 files for the choosen neuron
    for i = 1:length(dir)
        if(~isempty(strfind(content(i).name, 'msq1D.sa0')) && isempty(strfind(content(i).name, '.sa0.')))
            output(exp_index).events = fget_spk(content(i).name);
            output(exp_index).hdr = fget_hdr(content(i).name);
            output(exp_index).rate = length(output(exp_index).events)./total_time;
            exp_index = exp_index + 1;
        end
        if(~isempty(strfind(content(i).name, 'msq1d.sa0')) && isempty(strfind(content(i).name, '.sa0.')))
            output(exp_index).events = fget_spk(content(i).name);
            output(exp_index).hdr = fget_hdr(content(i).name);
            output(exp_index).rate = length(output(exp_index).events)./total_time;
            exp_index = exp_index + 1;
        end    
    end

    cd ../../..

catch
   disp('error'); 
end
end