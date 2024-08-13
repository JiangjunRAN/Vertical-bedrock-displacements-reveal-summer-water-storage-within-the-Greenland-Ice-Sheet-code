
# README


## 1. Introduction
This repository contains the core script for estimating the summer water storage from GNSS data.
The sample code is written in MATLAB. With the following file structure, you can run the demo directly.


## 2. File Structure
**Input**:&nbsp; .mat files for GNSS data<br>
&nbsp;<br>
**Src**<br>
　&nbsp; **--/inversion_Reg0E2013.m**：Main script to estimate the summer water storage.<br>
　&nbsp; **--/est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_Reg0ErMean_E_r_pavel.m**:&nbsp; Function used in inversion_Reg0E2013.m.<br>
　&nbsp;  **--/selected_time_series_to_given_time_interval.m**: &nbsp;Function used in inversion_Reg0E2013.m.<br>
 &nbsp;<br>
**Output**:&nbsp; Path to store the output.						 

## Corresponding authors
If you have any question about the codes , please contact Jiangjun RAN (ranjj@sustech.edu.cn).<br>
<br>
_Copyright (c) 2023 Jiangjun RAN. All rights reserved._ <br>
_Last update 2024.08.13_
