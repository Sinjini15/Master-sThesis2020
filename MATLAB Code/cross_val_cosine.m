clear all;
clc;

filename = "/Users/sinjinim/Documents/ASU Semesters/Thesis Sem/My code/val_set.csv"
peak_data = readtable(filename);
peak_data = removevars(peak_data,'Var1');
peak_data = table2array(peak_data);

start = find(peak_data>0,1);
stop = size(peak_data,1);
peak_data = peak_data(start:stop,:);
Fs = 500;
alpha= 0.05;

%%
%=== I: Scaling the data =====%

peak_data(:,1) = peak_data(:,1)/Fs;
peak_data(:,1) = peak_data(:,1) - min(peak_data(:,1));

training_data(:,1) = peak_data(:,1);
training_data(:,2) = peak_data(:,3);

%%
%=== II: Initializing variables =====%

loocv_data = training_data;
n= length(loocv_data(:,2));
h = linspace(0.01,35,2000);
%%
%=== III: LOOCV estimator =====%


tpts = loocv_data(:,1);
tpts= tpts';
rpts = loocv_data(:,2);
rpts= rpts';

for j = 1:length(tpts)
  t_pts(:,j) = tpts(setdiff(1:length(tpts),j));
  r_pts(:,j) = rpts(setdiff(1:length(rpts),j));
end

Wxseq = repmat(tpts,[n-1,1]);
Wyseq = repmat(rpts,[n-1,1]);

est2d_sum = zeros(1,length(h)); %to store the value of the sum of estimators for a specific h
avg = zeros(1, length(h)); % to store the final LLN average value for a specific h
%%

%---- Cosine Kernel

for i =1:length(h)
  
        W1 = (Wxseq - t_pts)/h(i);
        Wx = (pi/4).*(cos(pi/2 .* W1));
        W2 = (Wyseq - r_pts)/h(i);
        Wy = (pi/4).*(cos(pi/2 .* W2));
        est2d = Wx*Wy';
        est2d_sum (i) = sum(sum(est2d));
         
end

%%
%---Calculating the average 

for i=1:length(h)   
    avg(i)= 2*est2d_sum(i)/(n*(n-1)*h(i));
 end
%%
%=== IV: Calculation of integral/Convolution kernel =====%

pdf = zeros(1,length(h)); %to store the final values of the integral
int_pdf_sum = zeros(1,length(h)); %to store intermediate values

% Defining the matrices of data points used to calculate the kernel

Xi_tp = repmat(tpts,[n,1]);
Xi_rp = repmat(rpts,[n,1]);
Xj_tp = Xi_tp';
Xj_rp = Xi_rp';

%--- Cosine Kernel

for i = 1:length(h)
    tp_diff = abs((Xi_tp - Xj_tp)./h(i));
    rp_diff = abs((Xi_rp - Xj_rp)./h(i));
    tp_kern = ((pi/4).*cos(pi/2.*tp_diff)).^2;
    rp_kern = ((pi/4).*cos(pi/2.*rp_diff)).^2;
    int_pdf = tp_kern * rp_kern';
    int_pdf_sum(i) = sum(sum(int_pdf));
end

for i=1:length(h)
    
    pdf(i)=int_pdf_sum(i)/(n*n*h(i));
end

%%
%=== V: Calculation of estimated risk =====%

risk = pdf - avg;
[val,i] = min(risk);
h_final = h(i);


save h_val.mat h_final;

