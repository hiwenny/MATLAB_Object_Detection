function sample_test
clc; clf; clear all; close all;
[file_name, file_path] = uigetfile('*');

if file_name == 0
    return;
end

% im = imread([file_path, file_name]);
% % im = im(:,end:-1:1,:);
% mkdir('submission')
% filename = sprintf('%s.txt', file_name);
% z3351846_MTRN4230_ASST1([file_path, file_name], output_file_name, output_file_path);
% output in .txt
% cd ..
% % 
z3351846_MTRN4230_ASST1([file_path, file_name], file_name, '.\submission\IMG_0002.txt', '.auto_mark\student_submissions\z3351846_MTRN4230_ASST1\');