ANTERO-Gripper
MATLAB code for modeling, optimizing, and experimentally evaluating the ANTERO gripper and its 5-bar, 3-DOF (5B3D) finger mechanism.
This repository brings together:
a static finger model for contact-force estimation,
finger geometry and compliance optimization studies,
beam spring optimization scripts,
force-control analysis utilities, and
grasp-force validation scripts based on experimental data.
Repository structure
`5B3D-Finger-Model`
Baseline 5B3D finger static model and contact-force estimation code.
Key files:
`Main.m` — runs the main static-model sweep and saves results into the `Results/` directory.
`contactSolver.m`, `fastContactSolver.m`, `getJacobians.m`, `getParam.m`, `getSpring.m`, `updateParam.m` — supporting model functions.
`Results/` — generated outputs and follow-on plotting/analysis scripts.
`5B3D-Finger-Optimization`
Finger design-space exploration and optimization utilities.
Key files:
`Main.m` — runs a model sweep over finger states.
`Optimization.m` — performs design optimization over parameters such as link length and compliance-related terms.
`OptimizationAnalysis.m` — post-processes optimization outputs.
`FingerLinkOptimization.m` — helper script for optimization studies.
`Results/` — saved outputs and visualizations.
`5B3D-Spring-Optimization`
Beam spring modeling and stiffness optimization.
Key files:
`OptimizeSpringStiffness.m` — main entry point for spring optimization.
`solveBVP.m`, `bcfcn.m`, `initFcn.m` — boundary-value-problem setup.
`plotSpring.m`, `plotStress.m`, `plotStiffness.m`, `plotThickness.m` — visualization utilities.
`OptimalSpring.mat`, `OptimalSpring_Plus.mat` — saved optimization outputs.
`Force-Control-Testing`
Utilities for unpacking experimental datasets and analyzing force-control / frequency-response behavior.
Key files:
`unpackTeensy.m`
`unpackRokubi.m`
`interpData.m`
`FrequencyResponseAnalysis.m`
`Grasp-Force-Testing`
Scripts for comparing modeled grasp force against measured force-torque data.
Key files:
`Data/unpackTeensyData.m`
`Data/unpackRokubiData.m`
`interpData.m`
`GraspForceTestingAnalysis.m`
plotting utilities and exported figures
`ANTERO-Mid-Level-Control`
Placeholder folder for mid-level control content.
Requirements
This repository is written in MATLAB.
Depending on the scripts you run, you may need MATLAB functionality for:
symbolic math,
parallel computing, and
constrained optimization.
Several scripts also expect local experimental datasets to be available as `.csv`, `.log`, or `.mat` files in specific folder layouts.
Quick start
1. Clone the repository
```bash
git clone https://github.com/aambrose3/ANTERO-Gripper.git
cd ANTERO-Gripper
```
2. Open the repo in MATLAB
Open MATLAB and set the current folder to the repository root, or add the relevant subfolder to your MATLAB path.
3. Run the model or analysis you need
Static 5B3D finger model
```matlab
cd('5B3D-Finger-Model')
run('Main.m')
```
Finger optimization studies
```matlab
cd('5B3D-Finger-Optimization')
run('Main.m')              % model sweep
run('Optimization.m')      % design optimization
```
Spring optimization
```matlab
cd('5B3D-Spring-Optimization')
run('OptimizeSpringStiffness.m')
```
Grasp-force validation workflow
From `Grasp-Force-Testing`, run the scripts in this order:
```matlab
run('Data/unpackTeensyData.m')
run('Data/unpackRokubiData.m')
run('interpData.m')
run('GraspForceTestingAnalysis.m')
```
Expected workflow
5B3D model and optimization
Use the `5B3D-Finger-Model` and `5B3D-Finger-Optimization` folders when you want to:
sweep the finger through joint configurations,
estimate proximal/distal contact forces,
study feasible workspaces,
evaluate mechanical advantage, or
optimize design variables.
Spring design
Use `5B3D-Spring-Optimization` to generate and evaluate compliant spring designs and saved spring models.
Experimental validation
Use `Force-Control-Testing` and `Grasp-Force-Testing` to:
unpack Teensy and Rokubi datasets,
synchronize or interpolate time series,
analyze frequency response, and
compare modeled versus measured grasp force.
Important notes
Hard-coded paths and filenames
Many scripts use hard-coded relative paths and dataset names. Before running them, you will likely need to update file paths for your local machine and dataset organization.
Examples include paths such as:
`60x75\5_HZ\T1.csv`
`60x75\5_HZ\T1-FT.csv`
`Grasp_Force_Testing_Rokubi\90x90\T1.log`
`Data/Matched_Data_60x45.mat`
Results directories
Some scripts save outputs directly into `Results/` or into `.mat` files in the working directory. Running scripts from the intended folder helps avoid path issues.
Research code
This repository appears to be organized as a research codebase rather than a polished MATLAB package. Expect to inspect the top of each script for:
user settings,
sweep ranges,
optimization bounds,
dataset filenames, and
plotting/export behavior.
Suggested order for new users
If you are exploring the repo for the first time, a good path is:
Start with `5B3D-Finger-Model/Main.m` to understand the base model.
Move to `5B3D-Finger-Optimization/Optimization.m` for design studies.
Use `5B3D-Spring-Optimization/OptimizeSpringStiffness.m` for compliant spring design.
Run the experimental workflows in `Grasp-Force-Testing` and `Force-Control-Testing` once your data files are in place.
License
This project is licensed under the Apache License 2.0. See the `LICENSE` file for details.
Author
Repository author: Alexander B. Ambrose
---
If you use this code in academic work, consider citing the repository and any associated publications from the project.
