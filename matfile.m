% Use this a template to read files that are sequentially numbered. 

for i = 1:3

% Generate strings for file and image names

	  file_name = sprintf('inflammation-%02d.csv', i);
          img_name = sprintf('patient_data-%02.png', i);

          patient_data = csvread(file_name);
          ave_inflammation = mean(patient_data, 1);

          figure(i);

          subplot(2,2,1);
          plot(ave_inflammation);
          ylabel('average')

subplot(2,2,2);
plot(max(patient_data, [], 1));
ylabel('max')

subplot(2,2,3);
plot(min (patient_data, [], 1));
ylabel('min')

print('-dpng',img_name);
close();

end