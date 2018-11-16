warning off
clear
clc

load Data/Stimulus_Files/msq1D.mat
addpath('./MatlabFunctions/tview','./MatlabFunctions/fileload', './Data/Stimulus_Files','./Data/Spike_and_Log_Files', '.');

frames = length(msq1D);

files = dir ('./Data/Spike_and_Log_Files'); % Neuron Names
%%load newstruct.mat
load newstruct.mat

% These files have a problem and cannot be opened 
delete './Data/Spike_and_Log_Files/000620.a02/000620.a02lmsq1D.sa0'
delete './Data/Spike_and_Log_Files/000620.a02/000620.a02lmsq1D.log'

%% Section 2 - Question 3

j = 0;
disp('Neurons with spike count rate less than 2');
for i = 3 : length(files)
    j = j + 1;
    neurons(j).name = files(i).name;
    [neurons(j).data, neurons(j).freq] = Func_ReadData(files(i).name);
    neurons(j).avg_rate = mean([neurons(j).data.rate]);
    if(neurons(j).avg_rate < 2)
        disp(files(i).name);       % Unacceptable neurons
        neurons(j).acceptable = 0;
    else
        neurons(j).acceptable = 1;
    end
    fclose('all');
end

figure
histogram([neurons.avg_rate], 15);
title('Average Spike Count');


%% Section 2 - Question 5

cd('./Data/Spike_and_Log_Files/010718.B.c01');
tview('010718.B.c01atune.log');
cd ../../..

%% Section 3 - Question 1 and Question 2*

% The for loops are now removed from the code for the graders comfort!
% The sample neuron for our calculations
% Neuron Code = 000601.c05
% Neuron ID = 16 ;   neurons(16)
i = 16;

% Concating all neuron spike times in every experiment on the neuron
spike_times = [];
for j = 1 : length(neurons(i).data)
    a = neurons(i).data(j).events;
    spike_times = [spike_times, a'];
end

% Spike Triggered Stimuli Matrix:
neurons(i).stimuli = Func_StimuliExtraction (spike_times, msq1D, neurons(i).freq);

% Spike Triggered Average Method:
neurons(i).STA = mean(neurons(i).stimuli, 3); 

% Question 2*
% Performing TTest, Null Hyp. : Mean = 0
[~, p, ~] = ttest(permute(neurons(i).stimuli, [3, 1, 2]));
p = reshape(p , 16, 16);

figure
subplot(2,2,1);
imshow(neurons(i).STA, [-1, 1]);   % Rec Field with Normal Contrast
ylabel('Temporal'); xlabel('Spatial'); title([neurons(i).name,'   STA']);

subplot(2,2,2);
imshow(neurons(i).STA, [-0.1, 0.1]); % Rec Field with High Contrast
ylabel('Temporal'); xlabel('Spatial'); title('High Contrast');

subplot(2,2,[3, 4]);
imshow(1-p)                              %  T Test
ylabel('Temporal'); xlabel('Spatial'); title('T Test 1 - PValue');


%% Section 3 - Question 3


% Spike Triggered Stimuli's Projection on STA
s = size(neurons(i).stimuli);  
zeta = zeros(1, s(3));
for j = 1 : s(3)
   z = neurons(i).stimuli(:,:,j).*(neurons(i).STA);
   zeta(j) = sum(sum(z));
end

% Spike Triggered Stimuli Histogram
figure 
hold on
histogram(zeta, 'BinMethod', 'scott', 'Normalization', 'pdf');

% Random Spike Generation
L = length(zeta);
total_time = frames/neurons(i).freq;
random_spikes = 10000*total_time*rand(1, L);    
random_stimuli = Func_StimuliExtraction (random_spikes, msq1D, neurons(i).freq);

s = size(random_stimuli);
zeta_random = zeros(1, s(3));
for j = 1 : s(3)
   z1 = random_stimuli(:,:,j).*(neurons(i).STA);
   zeta_random(j) = sum(sum(z1));
end

% Random Spike Histogram
histogram(zeta_random,  'BinMethod', 'scott', 'Normalization', 'pdf')
legend('spike', 'control')
title(neurons(i).name)

% Fitting two Normal Distro.'s and Finding the Intersections 
x = linspace(min(zeta),max(zeta),10000);
g_real = pdf('Normal', x, mean(zeta), std(zeta));
g_random = pdf('Normal', x,  mean(zeta_random), std(zeta_random));
[~, index] = min(abs(g_random./g_real-1));
intersect = x(index);


y1 = get(gca,'ylim');
plot([intersect intersect], y1, 'LineWidth', 1.2, 'Color', 'yellow');
plot(x, g_real, 'LineWidth', 1.2, 'Color', 'blue')
plot(x, g_random,'LineWidth', 1.2,  'Color', 'red')

%% Section 3 - Question 4

% Performing T-Test on zeta and zeta_random
% Null Hyp: mean(zeta) = mean(zeta_random)

[~, p] = ttest2(zeta,zeta_random);
disp(['TTest P-Value: ', num2str(p)])

%% Section 3 - Question 5

% Percentage of correct predictions using the intesection point in S3-Q3
estimation = zeta > intersect;
neurons(i).STA_percent = (sum(estimation)/length(zeta))*100;
disp(['STA Correct Spike Predictions = ', num2str(neurons(i).STA_percent ), '%']);

%% Section 4 - Question 1


% The for loops are now removed from the code for the graders comfort!
% The sample neuron for our calculations
% Neuron Code = 000601.c05
% Neuron ID = 16 ;   neurons(16)
i = 16;


spike_times = [];
for j = 1 : length(neurons(i).data)
    a = neurons(i).data(j).events;
    spike_times = [spike_times, a'];
end

neurons(i).stimuli = Func_StimuliExtraction (spike_times, msq1D, neurons(i).freq);
s = size(neurons(i).stimuli);
stimuli = reshape(neurons(i).stimuli, 16*16, s(3));  

% The Spike Triggered Correlation Matrix
correltaion = corr(stimuli');
[eig_vec, eig_vals] = eig(correltaion, 'vector');


[neurons(i).eig_vals, ind] = sort(eig_vals);

% The most significant eigenvectors
neurons(i).vec1 = reshape(eig_vec(:, ind(256)), 16, 16);
neurons(i).vec2 = reshape(eig_vec(:, ind(255)), 16, 16);
neurons(i).vec3 = reshape(eig_vec(:, ind(254)), 16, 16);

figure
subplot(1,3,1);

imshow(neurons(i).vec1, [-0.2, 0.2]);
ylabel('Temporal'); xlabel('Spatial'); title([neurons(i).name, ' # 1']);

subplot(1,3,2);
imshow(neurons(i).vec2, [-0.2, 0.2]);
ylabel('Temporal'); xlabel('Spatial'); title([neurons(i).name, ' # 2']);

subplot(1,3,3);
imshow(neurons(i).vec3, [-0.2, 0.2]);
ylabel('Temporal'); xlabel('Spatial'); title([neurons(i).name, ' # 3']);


%% Section 4 - Question 2

L = length(neurons(i).stimuli);
total_time = frames/neurons(i).freq;

% Generating Random Spike Train and their corresponding eigenvalue and
% eigen vector:

random_eig_val = zeros(256, 20);
for j = 1 : 20
    random_spikes = 10000*total_time.*rand(1, L);                                      % The random spike train    
    random_stimuli = Func_StimuliExtraction (random_spikes, msq1D, neurons(i).freq);   
    s = size(random_stimuli);
    random_corr = corr(reshape(random_stimuli, s(3), 16 * 16));
    random_eig_val(:,j) = eig(random_corr, 'vector');
end
[random_eig_val, random_ind] = sort(random_eig_val);


random_eig_avg = mean(random_eig_val, 2);
random_eig_std = std(random_eig_val, 1, 2);


% Plotting the eigenvalues and confidence interval
figure
hold on
plot((random_eig_avg(226:256)) + 10.4*random_eig_std(226:256), '--');
plot((random_eig_avg(226:256)) - 10.4*random_eig_std(226:256), '--');
plot( neurons(i).eig_vals(226:256), 'Marker', '.')

legend('Mean + 10.4 STD', 'Mean - 10.4 STD', 'Real Eigenvalues');
title(neurons(i).name)
ylabel('eigenvalue')
xlabel('Rank')


%% Section 4 - Question 3
% Spike Triggered Stimuli's Projection on Eigenvectors
clear z zeta1 zeta2 zeta3 random_z random_zeta1 random_zeta2 random_zeta3

% Projection of spike triggered stimuli on eigenvecs
s = size(neurons(i).stimuli);  
zeta1 = zeros(1,s(3)); zeta2 = zeros(1,s(3)); zeta3 = zeros(1,s(3));
for j = 1 : s(3)
   z = neurons(i).stimuli(:,:,j) .* neurons(i).vec1;
   zeta1(j) = sum(sum(z));
   z = neurons(i).stimuli(:,:,j) .* neurons(i).vec2;
   zeta2(j) = sum(sum(z));
   z = neurons(i).stimuli(:,:,j) .* neurons(i).vec3;
   zeta3(j) = sum(sum(z));
end

% Projection of random stimuli on eigenvecs
random_spikes = 10000*total_time.*rand(1, 3*length(neurons(i).stimuli));                                      % The random spike train    
random_stimuli = Func_StimuliExtraction (random_spikes, msq1D, neurons(i).freq);
s = size(random_stimuli);
random_zeta1 = zeros(1,s(3)); random_zeta2 = zeros(1,s(3)); random_zeta3 = zeros(1,s(3));
for j = 1 : s(3)
   random_z = random_stimuli(:,:,j) .* neurons(i).vec1;
   random_zeta1(j) = sum(sum(random_z));
   random_z = random_stimuli(:,:,j) .* neurons(i).vec2;
   random_zeta2(j) = sum(sum(random_z));
   random_z = random_stimuli(:,:,j) .* neurons(i).vec3;
   random_zeta3(j) = sum(sum(random_z)); 
end

% Plotting the histograms
figure
subplot(2,1,1)
histogram(zeta2, 'Normalization' , 'pdf', 'BinMethod', 'scott');
hold on
histogram(random_zeta1, 'Normalization' , 'pdf');
legend('spike', 'control'); title([neurons(i).name, '  Eigenvector #1']); xlabel('\zeta_1');

subplot(2,1,2)
histogram(zeta1, 'Normalization' , 'pdf', 'BinMethod', 'scott');
hold on
histogram(random_zeta2, 'Normalization' , 'pdf');
legend('spike', 'control'); title([neurons(i).name, '  Eigenvector #2']); xlabel('\zeta_2');

figure
histogram2(zeta1, zeta2, 'Normalization' , 'pdf');
hold on
histogram2(random_zeta1,random_zeta2, 'Normalization' , 'pdf');
legend('spike', 'control'); title(neurons(i).name); xlabel('\zeta_1'); ylabel('\zeta_2');

%% Section 4 - Question 5

% 2D jointly gaussian distro.
a1 = -4:.2:4; a2 = -4:.2:4;
[X1,X2] = meshgrid(a1,a2);
mu = [mean(zeta1), mean(zeta2)];
Sigma = cov(zeta1, zeta2);
random_mu = [mean(random_zeta1), mean(random_zeta2)];
random_Sigma = cov(random_zeta1, random_zeta2);
F = mvnpdf([X1(:) X2(:)], mu ,Sigma);
F = reshape(F, length(a2), length(a1));
random_F = mvnpdf([X1(:) X2(:)],random_mu,random_Sigma);
random_F = reshape(random_F,length(a2),length(a1));

% 3D jointly gaussian distro.
a1 = -4:.2:4; a2 = -4:.2:4; a3 = -4:.2:4;
[X1,X2, X3] = meshgrid(a1,a2, a3);
mu3 = [mean(zeta1), mean(zeta2), mean(zeta3)];
ZETA3 = [zeta1; zeta2; zeta3]';
Sigma3 = cov(ZETA3);
random_mu3 = [mean(random_zeta1), mean(random_zeta2), mean(random_zeta3)];
random_Sigma3 = cov([random_zeta1; random_zeta2; random_zeta3]');
F3 = mvnpdf([X1(:) X2(:) X3(:)], mu3 ,Sigma3);
F3 = reshape(F3, length(a3), length(a2), length(a1));
random_F3 = mvnpdf([X1(:) X2(:) X3(:)],random_mu3,random_Sigma3);
random_F3 = reshape(random_F3,length(a3),length(a2),length(a1));

% Predicting Spikes by Finding Intersections of Three Normal Distros.
count3 = 0;
for k = 1 : length(zeta1)
    [~, index1] = min(abs(a1-zeta1(k)));
    [~, index2] = min(abs(a2-zeta2(k)));
    [~, index3] = min(abs(a3-zeta3(k)));
    if(F3(index1, index2, index3) > random_F3(index1, index2, index3))
        count3 = count3 + 1;
    end
end
% Predicting Spikes by Finding Intersections of Two Normal Distros.
count = 0;
for k = 1 : length(zeta1)
    [~, index1] = min(abs(a1-zeta1(k)));
    [~, index2] = min(abs(a2-zeta2(k)));
    if(F(index1, index2) > random_F(index1, index2))
        count = count + 1;
    end
end

if(count3 > count) 
    neurons(i).STC_percent = count3/length(zeta1);
    disp('3 Significant Eigenvectors');
    disp(['STC Correct Spike Predictions = ', num2str(100*count3/length(zeta1))]); 
end
if(count > count3) 
    neurons(i).STC_percent = count/length(zeta1);
    disp('2 Significant Eigenvectors');
    disp(['STC Correct Spike Predictions = ', num2str(100*count/length(zeta1))]); 
end


%% Section 5 - Question 1

for i = 1 : length(neurons)
    if ( neurons(i).acceptable)
        spike_times = [];
        for j = 1 : length(neurons(i).data)
            a = neurons(i).data(j).events;
            spike_times = [spike_times, a'];
        end

        neurons(i).stimuli = Func_StimuliExtraction (spike_times, msq1D, neurons(i).freq);
        s = size(neurons(i).stimuli);
        stimuli = reshape(neurons(i).stimuli, 16*16, s(3));  
        
       
        % The Spike Triggered Covariance Matrix
        covar = cov(stimuli');

        STAW = inv(covar)*reshape(neurons(i).STA, 256, 1);
        STAW = reshape(STAW, 16, 16);
         figure
         subplot(2,1,1)
         imshow(STAW, [-0.1, 0.1]);
         title([neurons(i).name, ' STA_w'])
         subplot(2,1,2)
         imshow(neurons(i).STA, [-0.1, 0.1]);
         title([neurons(i).name,' STA'])

        % Spike Triggered Stimuli's Projection on STA
        s = size(neurons(i).stimuli);  
        zeta = zeros(1, s(3));
        for j = 1 : s(3)
           z = neurons(i).stimuli(:,:,j).*(STAW);
           zeta(j) = sum(sum(z));
        end

        % Spike Triggered Stimuli Histogram
        figure 
        hold on
        histogram(zeta, 'BinMethod', 'scott', 'Normalization', 'pdf');

        % Random Spike Generation
        L = length(zeta);
        total_time = frames/neurons(i).freq;
        random_spikes = 10000*total_time*rand(1, L);    
        random_stimuli = Func_StimuliExtraction (random_spikes, msq1D, neurons(i).freq);

        s = size(random_stimuli);
        zeta_random = zeros(1, s(3));
        for j = 1 : s(3)
           z1 = random_stimuli(:,:,j).*(neurons(i).STA);
           zeta_random(j) = sum(sum(z1));
        end

        % Random Spike Histogram
        histogram(zeta_random,  'BinMethod', 'scott', 'Normalization', 'pdf')
        legend('spike', 'control')
        title(neurons(i).name)

        % Fitting two Normal Distro.'s and Finding the Intersections 
        x = linspace(min(zeta),max(zeta),10000);
        g_real = pdf('Normal', x, mean(zeta), std(zeta));
        g_random = pdf('Normal', x,  mean(zeta_random), std(zeta_random));
        [~, index] = min(abs(g_random./g_real-1));
        intersect = x(index);

        estimation = zeta > intersect;
        neurons(i).STAw_percent = (sum(estimation)/length(zeta))*100;
        disp(['STAw Correct Spike Predictions = ', num2str(neurons(i).STAw_percent ), '%']); 
        
    end
end
%%
close all
STA_P = [neurons.STA_percent];
STAW_P = [neurons.STAw_percent];

plot(STA_P, STAW_P, '.')
ylabel('STA_W Percent');
xlabel('STA Percent');
title('STA and STA_w comparison');

%% Section 5 - Question 2

% PRojection of 256 eigenvectors
i = 16;

s = size(neurons(i).stimuli);
stimuli = reshape(neurons(i).stimuli, 16*16, s(3));  

% The Spike Triggered Correlation Matrix
correltaion = corr(stimuli');
[eig_vec, eig_vals] = eig(correltaion, 'vector');

s = size(neurons(i).stimuli);
fullzeta = zeros(256, s(3));
for k = 1 : 256
    for j = 1 : s(3)
       z = neurons(i).stimuli(:,:,j) .* reshape(eig_vec(:,k), 16, 16);
       fullzeta(k,j) = sum(sum(z));
     end
end

fullrandzeta = zeros(256, s(3));
total_time = length(msq1D)/(neurons(i).freq);
random_spikes = 10000*total_time.*rand(1, 3*length(neurons(i).stimuli));                                      % The random spike train    
random_stimuli = Func_StimuliExtraction (random_spikes, msq1D, neurons(i).freq);
for k = 1 : 256
    for j = 1 : s(3)
       z = random_stimuli(:,:,j) .* reshape(eig_vec(:,k), 16, 16);
       fullrandzeta(k,j) = sum(sum(z));
    end
end

fullcov = cov(fullzeta');
fullrandcov = cov(fullzeta');
fullmean = mean(fullzeta');
fullrandmean =  mean(fullrandzeta');
c = 0;
for p = 1 : length(fullzeta)    
    RANDPDF =  mvnpdf((fullzeta(:,p))', fullrandmean, fullrandcov);
    PDF     =  mvnpdf((fullzeta(:,p))', fullmean, fullcov);
    if (PDF > RANDPDF)
        c = c + 1;
    end
end

disp(['256 Eigenvector Method Spike Prediction %: ', num2str(100*c/length(fullzeta))])

