% In this script, six scatter plots are generated for the following pairs
%  - E and R_{-1.9482}
%  - E and SCI_{-3.6624}
%  - E and SO_{-1.7285}
%  - ?H and R_{-1.2383}
%  - ?H and SCI_{-2.3554}
%  - ?H and SO_{-1.1165}
% ... with their respective regression lines.

close all; % Before drawing, close any figures already opened
clear;     % Clear all variables

% CONFIG: Line width and font fize for each curve in drawn figures
lineWidth = 2;
fontSize = 26;
% Save plots to images? Set to true if yes
saveToFile = false;
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
  "E", "?H" % Their labels
};

% 30 by 3 array, number of edges in each dudv partition: (2,2), (2,3) then (3,3)
d_f = [
  6 6 6 7  6 8  7 8 9 6  6  7  8  8  7 10 9  9  9  8  9  8 8 8  7 10  7  6  6  6
  0 4 8 6 12 8 10 8 6 8 16 14 12 12 14 8 10 10 10 12 10 12 8 8 10 12 10 20 12 16
  0 1 2 3  3 5  4 5 6 5  4  5  6  6  5 8  7  7  7  6  7  6 8 8  7  9 10  5 12 19
]';

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

for edn = 1:numData % edn = experimental data number | 1=E, 2=?H
  for fnn = 1:numIndices % fnn = function number | 1=R_a, 2=SCI_a, 3=SO_a
    ccFn = @(alpha) corrcoef( % Gets corrcoef between benzenoid prop and index
      getIndexFns{fnn}(alpha)(!isnan(expData{1,edn})),
      expData{1,edn}(!isnan(expData{1,edn}))
    )(1,2);

    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, -5, 5, 1e-15));

    % Regression line i.e., y = mx + b;
    model = polyfit(getIndexFns{fnn}(peakAlpha), expData{1,edn},1);
    m = model(1); b = model(2);
    x = [getIndexFns{fnn}(peakAlpha)(1) getIndexFns{fnn}(peakAlpha)(end)];
    y = m*x + b;

    % Scatter plot
    this_figure = figure(3*(edn-1)+fnn); hold on;
    regLine = plot(x, y, '-', 'LineWidth', lineWidth);
    points = plot(getIndexFns{fnn}(peakAlpha), expData{1,edn}, '*', 'MarkerSize', 8, 'LineWidth', lineWidth/2);
    bestIndexLabel = sprintf("%s_{?%s}", indexName{fnn}, as_4_dp_str(abs(peakAlpha)));
    pointsLabel = sprintf("%s and %s", bestIndexLabel, expData{2,edn});

    % Label the scatter plot
    title(sprintf('between %s and %s', expData{2,edn}, bestIndexLabel));
    xlabel(bestIndexLabel);
    ylabel(sprintf('%s', expData{2,edn}));
    xticklabels(strrep(xticklabels,'-','?'));
    yticklabels(strrep(yticklabels,'-','?'));
    leg = legend("Regression Line", "Actual Data");
    set(leg, 'location', "southeast");

    % Change the font size to size set in the start of the script
    set(findall(this_figure,'-property','FontSize'),'FontSize', fontSize)

    drawnow;
  end
end

if saveToFile
  % Also enforce more suitable axes for figures 2 and 5
  saveas(figure(1), "03_scatter_E_R.png");
  figure(2);
  axis([0.035 0.1125 250 600]);
  saveas(figure(2), "03_scatter_E_SCI.png");
  saveas(figure(3), "03_scatter_E_SO.png");
  saveas(figure(4), "03_scatter_DH_R.png");
  figure(5);
  axis([0.2 0.9625 50 400]);
  saveas(figure(5), "03_scatter_DH_SCI.png");
  saveas(figure(6), "03_scatter_DH_SO.png");
end