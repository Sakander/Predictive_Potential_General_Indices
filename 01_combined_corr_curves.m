% In this script, six values are closely approximated via golden section search
%  - α value for which correlation coefficient ρ is strongest between E and R_a
%  - α value for which ρ is strongest between E and SCI_α
%  - α value for which ρ is strongest between E and SO_α
%  - α value for which ρ is strongest between ΔH and R_α
%  - α value for which ρ is strongest between ΔH and SCI_α
%  - α value for which ρ is strongest between ΔH and SO_α
% Additionally, curves for ρ against α near these values are plotted in 4 figs

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 2;
fontSize = 16;
saveToFile = false; % Set to true to auto-save plots
% Note: Dimensions of resulting images are scaled according to each window size.
%       To control the dimensions, after running the whole script, resize each
%       ... figure window and then run only the saveas functions
%       ... manually, by selection, at the end of this script

% Utility: Below is a function to round-off to 4 decimal places | returns string
%          Need to use this function as round(X,4,Type) does not exist in Octave
%          ... and sprintf("%.04f",X) does not round properly for some numbers.
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));

% Cell containing Entropy and Heat Capacity of 30 lower benzenoids
expData = {reshape([ % Entropy
    269.722 334.155 389.475 395.882 444.724 447.437
    457.958 455.839 450.418 399.491 499.831 513.857
    508.537 507.395 506.076 512.523 500.734 510.307
    509.210 513.879 511.770 509.611 461.545 463.738
    468.712 555.409 472.295 554.784 468.796 551.708
  ]', 30, 1), reshape([ % Heat Capacity
    83.019 133.325 184.194 183.654 235.165 233.497
   234.568 234.638 233.558 200.815 286.182 285.056
   284.037 284.088 285.148 284.595 284.870 284.503
   284.785 284.740 284.233 284.552 251.175 250.568
   251.973 336.098 267.543 337.204 285.041 368.518
  ]', 30, 1);
  "E", "ΔH" % Their labels
};

% 30 by 3 array, number of edges in each dudv partition: (2,2), (2,3) then (3,3)
d_f = [
  6 6 6 7  6 8  7 8 9 6  6  7  8  8  7 10 9  9  9  8  9  8 8 8  7 10  7  6  6  6
  0 4 8 6 12 8 10 8 6 8 16 14 12 12 14 8 10 10 10 12 10 12 8 8 10 12 10 20 12 16
  0 1 2 3  3 5  4 5 6 5  4  5  6  6  5 8  7  7  7  6  7  6 8 8  7  9 10  5 12 19
]'; % Used for computing indices based on edge endpoint degree partitions

% Cell containing the three index-computing functions
% Usage: getIndexFns{n}(a) | n=1:R_a, n=2:SCI_a, n=3:SO_a, a = alpha
getIndexFns = { % obtains a 30 by 1 matrix containing indices of the benzenoids
  @(a) (sum(d_f.*[4,6,9].^a,2)); % General Randic index
  @(a) (sum(d_f.*[4,5,6].^a,2)); % General SCI
  @(a) (sum(d_f.*[8,13,18].^a,2)); % General Sombor index
}';
% Cell containing their respective labels
indexName = {"R" "SCI" "SO"};

% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % two
numIndices = size(getIndexFns,2); % three
numCurves = numData*numIndices;   % six

% All x in visible ranges (both plots - near and far)
% E's xnear ∈ [-4,0.25]; ΔH's xnear ∈ [-3,0.25]
xnear = linspace(-1,0,800) .* [4.25; 3.25] + 0.25;
xfar = linspace(-40,40,800); % xfar is same for both E and ΔH

% Do the same procedure for each experimental data i.e., E and ΔH
for ii = 1:numData
  % This figure is for the zoomed-in plot
  figure(numData+ii); hold on;

  % Indicate inverval for which R_a is the better indicator (before corr-curves)

  % WARNING: these xmeet1-ymeet2 values are hardcoded, computed separately
  %          ... These are coordinates where ρ-α of {ΔH or E}-R_α
  %          ... intersects ρ-α of {ΔH or E}-SCI_α
  xmeet1 = [
    -2.57876284788673; % for E
    -1.65556221240371  % for ΔH
  ](ii); % of these 2, either the first value or second is used, depending on ii
  ymeet1 = [
    0.977567833002574; % for E
    0.994298524387525  % for ΔH
  ](ii);
  xmeet2 = 0;
  ymeet2 = [
    0.919815128725857; % for E
    0.974989505892947  % for ΔH
  ](ii);

  % Plot the blue dashed box (before drawing the curves so it appear beneath)
  plot([xmeet1 xmeet1 0 0], [0 1 1 0], '--b', 'LineWidth', lineWidth);

  yend = 0; % <-- to be assigned some value later for adjusting visible range

  % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    % Function to get corrcoef ρ between E/ΔH (depending on ii) with specified α
    %                                and R_a/SCI_a/SO_a (depending on n)
    ccFn = @(alpha) corrcoef(
      getIndexFns{n}(alpha)(!isnan(expData{1,ii})),
      expData{1,ii}(!isnan(expData{1,ii}))
    )(1,2);

    % generate corresponding y values
    ynear = arrayfun(ccFn, xnear(ii,:));
    yfar = arrayfun(ccFn, xfar);

    % Compute peak values via. golden section search, and display in console
    disp(sprintf("%s against general_%s", expData{2,ii}, indexName{n}));
    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, xnear(ii,1), xnear(ii,end), 1e-15))
    peakCorrCoeff = ccFn(peakAlpha)

    % Generate curve label                  [E,ΔH]         [R,SCI,SO]
    curveLabels{n} = sprintf("%s and %s_a", expData{2,ii}, indexName{n});

    figure(ii); % One zoomed-out plot for each expData
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curvesFar(n) = plot(xfar, yfar, '-', 'LineWidth', lineWidth);
    drawnow;

    figure(numData+ii); % Each expData's set of curves has a zoomed-in plot
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curves(n) = plot(xnear(ii,:), ynear, '-', 'LineWidth', lineWidth);

    % Show the peak in the plot: draw indicator lines & display coordinates
    plot([peakAlpha peakAlpha xnear(ii,1)], [0 peakCorrCoeff peakCorrCoeff],
         '--k', 'LineWidth', lineWidth/2); % Black dashed indicator lines
    text(peakAlpha, peakCorrCoeff,
        {'', sprintf("(−%s, %s)", as_4_dp_str(abs(peakAlpha)), as_4_dp_str(peakCorrCoeff))},
        'VerticalAlignment', 'bottom');
        % Negative sign entered manually here to bypass the default
        % ... usage of "hypen-minus" instead of "minus" (− vs -)

    yend = max(yend, ynear(end)); % y value to be used as visible y lower bound
  end

  % Mark and write on the plot the limits of alpha where R_a is better than SCI_a
  plot([xmeet1 xmeet2], [ymeet1 ymeet2], '*b',
       'MarkerSize', 16, 'LineWidth', lineWidth/1.5); % Mark with blue asterisks
  text(xmeet1, ymeet1, {'', sprintf(" (−%s, %s)", as_4_dp_str(abs(xmeet1)), as_4_dp_str(ymeet1))},
       'VerticalAlignment', 'top', 'Color', [0, 0, 0.8]); % Write blue text
  text(xmeet2, ymeet2, {'', sprintf("(0, %s) ", as_4_dp_str(ymeet2))},
       'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', [0, 0, 0.8]);

  % Label this expData's zoomed-in plot
  xlabel('α');
  ylabel('ρ');
  leg = legend(curves, curveLabels); % curves contains all drawn "xnear" curves
  set(leg, 'location', "southwest");
  axis([xnear(ii,1) xnear(ii,end) yend 1.001]); % Enforce figure's visible range
  drawnow;

  % Label the zoomed-out plot
  figure(ii);
  xlabel('α');
  ylabel('ρ');
  leg = legend(curvesFar, curveLabels); % curvesFar contains all drawn "xfar" curves
  if ii==2
    set(leg, 'location', "southeast");
  end

  hold off;
end

for ii = 1:4
  % Replace hyphens with minuses on negative axes
  figure(ii);
  xticklabels(strrep(xticklabels,'-','−'));
  yticklabels(strrep(yticklabels,'-','−'));

  % Set all fontsizes to size specified early in the script
  set(findall(figure(ii),'-property','FontSize'),'FontSize', fontSize)
end

if saveToFile
  saveas(figure(1), "01_comb_ccurves_E_indices_FAR.png");
  saveas(figure(2), "01_comb_ccurves_DH_indices_FAR.png");
  saveas(figure(3), "01_comb_ccurves_E_indices_NEAR.png");
  saveas(figure(4), "01_comb_ccurves_DH_indices_NEAR.png");
end
