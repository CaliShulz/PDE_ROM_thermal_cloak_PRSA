%% run_all.m 

clear all
close all
clc

% Run all the scripts needed to generate the figures

mesh_names = dir('*.mat') ;

for tt = 1:length(mesh_names)
   
   close all
   clc
   
   mesh_name = erase(mesh_names(tt).name,'.mat');
   
%    assemblingStart = tic();
%    assemble_problem(mesh_name);
%    assembling_time(tt) = toc(assemblingStart);
%    
%    SSreductionStart = tic();
%    SS_reduction(mesh_name);
%    SS_ocp_time(tt) = toc(SSreductionStart);
%    
%    trainingStart = tic();
%    T_POD_training(mesh_name);
%    T_POD_training_time(tt) = toc(trainingStart);
%    
%    T_reductionStart = tic();
%    T_reduction(mesh_name);
%    T_reduction_time(tt) = toc(T_reductionStart);
%    
%    T_reduction_testsStart = tic();
%    T_reduction_tests(mesh_name);
%    T_reduction_tests_time(tt) = toc(T_reduction_testsStart);
   
   T_plotsStart = tic();
   T_plots(mesh_name);
   T_plots_time(tt) = toc(T_plotsStart);
   
    
end

% sslash = path_setup() ; % setup path 
% file_id = fopen(strcat('archive_sim',sslash,'report','.txt'),'w')
% fprintf(file_id,'################################### \n');
% 
% for tt = 1:length(mesh_names)
% 
% fprintf(file_id,'\n');
% fprintf(file_id,'CASE: %s \n', erase(mesh_names(tt).name,'.mat'));
% fprintf(file_id,'Assembling Time : %d \n', assembling_time(tt));
% fprintf(file_id,'Steady State OCP Time : %d \n', SS_ocp_time(tt));
% fprintf(file_id,'Training Time : %d \n', T_POD_training_time(tt));
% fprintf(file_id,'Reduction Time : %d \n', T_reduction_time(tt));
% fprintf(file_id,'Reduction Tests Time : %d \n', T_reduction_tests_time(tt));
% fprintf(file_id,'Plotting Time: %d  \n' , T_plots_time(tt));
% 
% end

