% In this script, six plots are generated which shows α intervals for
%  - good correlation coefficient ρ between E and R_α
%  - good correlation coefficient ρ between E and SCI_α
%  - good correlation coefficient ρ between E and SO_α
%  - good correlation coefficient ρ between ΔH and R_α
%  - good correlation coefficient ρ between ΔH and SCI_α
%  - good correlation coefficient ρ between ΔH and SO_α
% The limits of the intervals are also printed to console

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 3;
fontSize = 24;
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

% Boundaries for visible intervals, for each index-property pair
%             R_a   SCI_a   SO_a
xstart = [   -3.2    -5.5   -2.2   % E
             -2.4    -4.4   -1.9]; % ΔH
xend = [     -0.5    -1.5   -1.2   % E
             0.25    0.25   -0.1]; % ΔH
ystart = [0.99301 0.98901  0.983;  % E
           1.0001  1.0001 0.9981]; % ΔH
yend = [     0.95    0.96  0.975;  % E
             0.97   0.975   0.98]; % ΔH

% Exact rho value considered good for
%            E     ΔH
a_goodrho = [0.98; 0.99];

% Colors (different shades of cyan and green)
colShaded = {[0.85, 1, 1]; [0.85, 1, 0.85]};
colIndicatorVert = {[0.2, 0.55, 0.55]; [0, 0.5, 0]};
colIndicatorHorz = {[0.35, 0.75, 0.75]; [0.45, 0.7, 0.45]};
colCurve = {[0, 0.75, 0.75]; [0, 0.5, 0]};

% Do the same procedure for each experimental data i.e., E and ΔH
for ii = 1:numData
  % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    % Function to get corrcoef ρ between E/ΔH (depending on ii) with specified α
    %                                and R_a/SCI_a/SO_a (depending on n)
    ccFn = @(alpha) corrcoef(
      getIndexFns{n}(alpha)(!isnan(expData{1,ii})),
      expData{1,ii}(!isnan(expData{1,ii}))
    )(1,2);

    this_fig = figure((ii-1)*numIndices+n); % Basically just figures (1) to (6)
    hold on;

    % Get Interval limits and print them to console

    % Get alpha for highest rho first
    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, -4, 0, 1e-15));

    % Prepare a funcion to calc |rho(a)-goodrho|
    ccFn_good = @(a)(-abs(ccFn(a)-a_goodrho(ii)));

    % and func to get the limit, i.e., value of alpha where rho is 0.98 or 0.99
    getLimitFromInterval = @(lb, ub) mean(
      GoldenSectionSearch_Maximum(ccFn_good, lb, ub, 1e-15));

    a_lb = getLimitFromInterval(peakAlpha-3, peakAlpha); % Search to the left
    a_ub = getLimitFromInterval(peakAlpha, peakAlpha+3); % Search to the right

    % Write the intervals in console
    disp(sprintf("ρ(%s,%s_α) ≥ %.02f when α ∈ [%.08f, %.08f]",
         expData{2,ii}, indexName{n}, a_goodrho(ii), a_lb, a_ub));

    % Plot the actual curve, but exclude good alpha range (<- separately drawn)
    % generate x values
    x = [linspace(xstart(ii,n),a_lb,300), linspace(a_ub,xend(ii,n),300)];
    % generate the corresponding y values
    y = arrayfun(ccFn, x);
    plot(x, y, '-', 'LineWidth', lineWidth);

    % Shade the area inside the good alpha interval
    u0 = a_lb;         u_width = a_ub-a_lb;
    v0 = ystart(ii,n); v_height = yend(ii,n) - ystart(ii,n);
    rectangle('Position', [u0, v0, u_width, v_height], 'FaceColor', colShaded{ii}, 'LineStyle', 'none');

    % Draw the indicator lines for the good alpha interval
    % Vertical dashed lines:
    plot([a_lb a_lb], [ystart(ii,n) yend(ii,n)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});
    plot([a_ub a_ub], [ystart(ii,n) yend(ii,n)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});
    % Horizontal black dashed line, horizontal colored dashed line:
    plot([xstart(ii,n) a_lb], [a_goodrho(ii) a_goodrho(ii)],
         '--k', 'LineWidth', lineWidth/1.75);
    plot([a_lb a_ub], [a_goodrho(ii) a_goodrho(ii)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorHorz{ii});

    % Write on the plot the interval limits
    text(a_lb, yend(ii,n), {'', sprintf("  α=−%s", as_4_dp_str(abs(a_lb)))}, 'VerticalAlignment', 'top', 'Rotation', 90);
    text(a_ub, yend(ii,n), {'', sprintf("  α=−%s", as_4_dp_str(abs(a_ub)))}, 'VerticalAlignment', 'top', 'Rotation', 90);

    % Finally, plot the colored curve within the good alpha range
    x_in = linspace(a_lb,a_ub,400);
    y_in = arrayfun(ccFn, x_in);
    plot(x_in, y_in, '-', 'LineWidth', lineWidth, 'Color', colCurve{ii});
    % Also highlight the good interval on the x-axis
    plot([a_lb, a_ub], [yend(ii,n), yend(ii,n)], '-',
         'LineWidth', lineWidth, 'Color', colCurve{ii});

    % Label the plot                     [E,ΔH]         [R,SCI,SO]
    title(sprintf('between %s and %s_a', expData{2,ii}, indexName{n}));
    xlabel('α');
    ylabel('ρ');
    drawnow;
    axis([xstart(ii,n) xend(ii,n) yend(ii,n) ystart(ii,n)]);

    % Replace all hypens with minuses
    xticklabels(strrep(xticklabels,'-','−'));
    yticklabels(strrep(yticklabels,'-','−'));
    % Set all fontsizes to size specified early in the script
    set(findall(this_fig,'-property','FontSize'),'FontSize', fontSize);
    drawnow;

    hold off;

  end
end

if saveToFile
  saveas(figure(1), "02_good_a_intervals_E_R_a.png");
  saveas(figure(2), "02_good_a_intervals_E_SCI_a.png");
  saveas(figure(3), "02_good_a_intervals_E_SO_a.png");
  saveas(figure(4), "02_good_a_intervals_DH_R_a.png");
  saveas(figure(5), "02_good_a_intervals_DH_SCI_a.png");
  saveas(figure(6), "02_good_a_intervals_DH_SO_a.png");
end
