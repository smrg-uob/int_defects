#int_defects.sh
#Harry Coules 2017
#	- Modified 03/2019: Added automatic write/deletion of wrapper functions.
#
#DESCRIPTION
#This Linux shell script is used for creating and executing a set of finite element analyses of structures containing crack-like defects.It boots Matlab and runs two scripts: int_defects_write_inp_wrapper.m creates a set of model input files using Abaqus/CAE, and int_defects_run_wrapper.m runs the models using Abaqus/Standard.
#
#INSTRUCTIONS
#1. Ensure that the int_defects toolbox is installed on the MATLAB installation which will be used, and the Abaqus package is available.
#2. Before running, edit this script so that the lines for creating "wrapper functions" contain correct parameters for the model set and execution conditions that you want to use.
#3. The script should be placed in a directory containing a master .py file for the model type which is being used and a Matlab .mat file containing a structure which defines the range of model parameters (paramRangeStruct).
#4. Run: chmod +x int_defects.sh	#Set execute permission
#		 bash int_defects.sh		#Execute int_defects.sh

#Create wrapper functions - THESE LINES WILL NEED TO BE EDITED DEPENDING ON YOUR APPLICATION
echo "int_defects_write_inp_parametric('paramRangeStruct_NASA_P1_2.mat','IntCrackJob1_6.0.4_single_master.py','abaqus',600);" >> int_defects_write_inp_wrapper.m
echo "load('abaqus_param_files/workspace_dump.mat'); int_defects_run_parametric_parallel(paramLogStruct,filenamesStruct,'abaqus',44000,1,4,false,true);" >> int_defects_run_wrapper.m

#Create and run models
nohup matlab -nodisplay <int_defects_write_inp_wrapper.m >matlab_cmd_log1.txt
nohup matlab -nodisplay <int_defects_run_wrapper.m >matlab_cmd_log2.txt

#Delete wrapper functions
rm int_defects_write_inp_wrapper.m
rm int_defects_run_wrapper.m