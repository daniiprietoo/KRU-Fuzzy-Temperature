%% 1. Initialize Fuzzy System
fis = mamfis('Name', 'Real_Feel_System');

% --- Configuration ---
% Ranges
t_range = [-10 40];   % Temperature (Celsius)
h_range = [0 100];    % Humidity (%)
out_range = [-15 60]; % Apparent Temp (Humidex)

% Granularity (Standard Precision)
num_t_mfs = 5;    % "VeryCold" ... "Hot"
num_h_mfs = 5;    % "Dry" ... "Wet"
num_out_mfs = 9;  % Higher resolution for output

%% 2. Generate Inputs (Gaussian MFs)
% Temperature Input
t_step = (t_range(2) - t_range(1)) / (num_t_mfs - 1);
t_sigma = t_step / 2.355;
fis = addInput(fis, t_range, 'Name', 'Temperature');

for i = 1:num_t_mfs
    center = t_range(1) + (i-1)*t_step;
    fis = addMF(fis, 'Temperature', 'gaussmf', [t_sigma center], 'Name', "T_Level"+i);
end

% Humidity Input
h_step = (h_range(2) - h_range(1)) / (num_h_mfs - 1);
h_sigma = h_step / 2.355;
fis = addInput(fis, h_range, 'Name', 'Humidity');

for i = 1:num_h_mfs
    center = h_range(1) + (i-1)*h_step;
    fis = addMF(fis, 'Humidity', 'gaussmf', [h_sigma center], 'Name', "H_Level"+i);
end

%% 3. Generate Output (Apparent Temp)
out_step = (out_range(2) - out_range(1)) / (num_out_mfs - 1);
out_sigma = out_step / 2.355;
fis = addOutput(fis, out_range, 'Name', 'ApparentTemp');

% Store centers for mapping
out_centers = zeros(1, num_out_mfs); 
for i = 1:num_out_mfs
    center = out_range(1) + (i-1)*out_step;
    out_centers(i) = center;
    fis = addMF(fis, 'ApparentTemp', 'gaussmf', [out_sigma center], 'Name', "Out_Level"+i);
end

%% 4. Rule Generation with "Desert Safety" Logic
ruleList = [];

% Loop through every combination of Temp and Humidity
for t_idx = 1:num_t_mfs
    T_val = t_range(1) + (t_idx-1)*t_step;
    
    for h_idx = 1:num_h_mfs
        H_val = h_range(1) + (h_idx-1)*h_step;
        
        % --- STEP A: Calculate Vapor Pressure (e) ---
        % Using the precise exponential formula
        exponent_term = (17.27 * T_val) / (237.7 + T_val);
        e_val = (H_val / 100) * 6.105 * exp(exponent_term);
        
        % --- STEP B: Calculate Humidex ---
        % Formula: T + (5/9)*(e - 10)
        raw_humidex = T_val + (5/9) * (e_val - 10);
        
        % --- STEP C: The "Desert Safety" Constraint ---
        % PROBLEM: In dry heat (e < 10), raw_humidex subtracts heat.
        % FIX: Ensure Feels-Like never drops more than 2 degrees below Actual.
        % This forces 40C/0% Hum to be at least 38C (not 31C).
        final_AT = max(raw_humidex, T_val - 2); 
        
        % --- STEP D: Map to Nearest Output MF ---
        [~, out_idx] = min(abs(out_centers - final_AT));
        
        % Add Rule: [Temp_Idx, Hum_Idx, Output_Idx, Weight, Operator]
        ruleList = [ruleList; t_idx h_idx out_idx 1 1];
    end
end

% Add the rules to the system
fis = addRule(fis, ruleList);

%% 5. Visualization and Validation
% Plot the surface to verify the "Dip" is gone
figure;
gensurf(fis);
title('Corrected Apparent Temperature (No Artificial Cooling)');
xlabel('Temperature');
ylabel('Humidity');
zlabel('Feels Like Temp');

% Open the Designer App to inspect manually
fuzzyLogicDesigner(fis)