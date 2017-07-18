function printMetrics(metrics, perJointMetrics, dispHeader, dispMetrics)
% print metrics
% 
% ...
%

metricsInfo.names.long = {'Recall','Precision','False Alarm Rate', ...
    'GT Tracks','Mostly Tracked','Partially Tracked','Mostly Lost', ...
    'False Positives', 'False Negatives', 'ID Switches', 'Fragmentations', ...
    'MOTA','MOTP', 'Head', 'Shoulders', 'Elbows', 'Wrists', 'Hips', 'Knees', 'Ankles'};

metricsInfo.names.short = {'Rcll','Prcn','FAR', ...
    'GT','MT','PT','ML', ...
    'FP', 'FN', 'IDs', 'FM', ...
    'MOTA','MOTP', 'Head', 'Sho ', 'Elb ', 'Wri ', 'Hip ', 'Knee', 'Ank '};

selected = {'Rcll','Prcn', ...
    'GT','MT','PT','ML', ...
    'FP', 'FN', 'IDs', 'FM', ...
    'MOTA','MOTP', 'Head', 'Sho ', 'Elb ', 'Wri ', 'Hip ', 'Knee', 'Ank '};

metricsInfo.widths.long = [6 9 16 9 14 17 11 15 15 11 14 6 5 8 4 9 6 6 4 5 6];
metricsInfo.widths.short = [5 5 5 3 3 3 3 5 5 4 4 5 6 5 5 5 5 5 5 5 5];

metricsInfo.format.long = {'.1f','.1f','.2f', ...
    'i','i','i','i', ...
    'i','i','i','i', ...
    '.1f','.1f','.1f', ...
    '.1f','.1f','.1f','.1f','.1f','.1f','.1f'};

metricsInfo.format.short=metricsInfo.format.long;    


namesToDisplay=metricsInfo.names.long;
widthsToDisplay=metricsInfo.widths.long;
formatToDisplay=metricsInfo.format.long;

namesToDisplay=metricsInfo.names.short;
widthsToDisplay=metricsInfo.widths.short;
formatToDisplay=metricsInfo.format.short;

if nargin<3, dispHeader=1; end

% padChar={' ',' ','|',' ',' ',' ','|',' ',' ',' ','| ',' ','|',' ', ' ',' ',' ',' ',' ',' ',' ',' '};
padChar={' & ',' & ',' & ',' & ',' & ',' & ',' & ',' & ',' & ',' & ',' & ',' & ',' & ',' & ', ' & ',' & ',' & ',' & ',' & ',' & ',' & ',' & '};


numMetrics = length(metrics) + length(perJointMetrics);

if nargin<4
    dispMetrics=1:numMetrics;
end

if dispHeader
    for m=dispMetrics
        if(m <= length(metrics))
            printString=sprintf('fprintf(''%%%is%s'',char(namesToDisplay(m)))',widthsToDisplay(m),char(padChar(m)));
            eval(printString)
        else
            printString=sprintf('fprintf(''%%%is%s'',char(namesToDisplay(m)))',2*widthsToDisplay(m),char(padChar(m)));
            eval(printString)
        end
    end
    fprintf('\n');
end

for m=dispMetrics
    if(m <= length(metrics))
        printString=sprintf('fprintf(''%%%i%s%s'',metrics(m))',widthsToDisplay(m),char(formatToDisplay(m)),char(padChar(m)));
        eval(printString)
    else
        i = m-length(metrics);
        printString=sprintf('fprintf(''%%%i%s/%%%s%s'',perJointMetrics(1,i),perJointMetrics(2,i))',widthsToDisplay(m),char(formatToDisplay(m)),char(formatToDisplay(m)),char(padChar(m)));
        eval(printString)
    end
end



% if standard, new line
if nargin<5
    fprintf('\n');
end