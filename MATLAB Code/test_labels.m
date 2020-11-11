
%== Load test dataset to extract the test set annotations ==%
filename = "/Users/sinjinim/Documents/ASU Semesters/Thesis Sem/My code/test_set.csv"
test_data = readtable(filename);
test_data = removevars(test_data,'Var1');
test_data = table2array(test_data);
start = find(test_data>0,1);
stop = size(test_data,1);
test_data = test_data(start:stop,:);
Fs = 500;
alpha= 0.05;


test_ann = test_data(:,2);

%%
%== Load the annotation file for infant==%
event_info = load('event_gen.mat');
ann = event_info.ann;
n = length(ann);
testn = length(test_ann);

%=== Check to see which of the test annotations are unhealthy ===%
% testing for if any of the test set annotations lie close to the ann file

for i=1:testn
    test_val = test_ann(i);
    for j=1:n
        ann_val = ann(j);
        u = rand;
        if test_val > ann_val && test_val < ann_val+(ceil(u*1500))
            unhealthy_ann(i,j) = 1;
        else
            unhealthy_ann(i,j)=0;
        end
    end
end

[row,col] = find(unhealthy_ann);

%%
%=== Generate which of the points in the test set are unhealthy ===%
row = sort(row);

save brad_labels.mat row;

        
