# Code for Hydraulic Optimization of Canal Regulation in Verziano Waste Pretreatment Plant

This repository contains the MATLAB scripts associated with a hydraulic study conducted on a canal conveying sewage to a waste pretreatment plant in Verziano, Brescia. The objective of the study was to design the regulation of devices that ensure a constant flow rate to the pretreatment plant, especially during periods of variable and high incoming discharge, such as during rainfall.

## Project Overview

The canal system at the Verziano pretreatment plant has several discharge regulation mechanisms, including a lateral weir and a sluice gate, that work together to manage varying incoming discharge and ensure the proper operation of the plant. When the incoming discharge exceeds the capacity of the plant, part of the sewage is diverted to a surface water body.

## MATLAB Scripts

This repository contains two MATLAB scripts used in the hydraulic study:

1. **Verziano_1: Lateral Weir Crest Height Calculation**  
   This script accepts one or more incoming discharge values as input and computes the amount by which the lateral gates need to be lowered. The output is the required depth to lower the gates for each input discharge, with results for one up to four gates in operation. If multiple gates are in operation, all gates are assumed to be lowered by the same amount.

2. **Verziano_2: Sluice Gate Opening Calculation**  
   This script determines the opening of the sluice gate for a fixed incoming discharge (12 m³/s) and a specific lateral weir crest height. The crest height is computed from the results of the first part of the study.

## How to Use

1. **Clone this repository**:  
   `git clone https://github.com/yourusername/repository-name.git`

2. **Run the Scripts**:  
   - Script 1: `lateral_weir_crest_height.m`
     - Input: A single or multiple discharge values.
     - Output: The required lowering of the gates for each discharge value.
   - Script 2: `sluice_gate_opening.m`
     - Input: The calculated lateral weir crest height and the incoming discharge of 12 m³/s.
     - Output: The corresponding sluice gate opening.

3. **Modify the scripts as needed** to adapt the calculations to different canal geometries or discharge conditions.

## Important Notes

- The scripts assume a fixed sluice gate height for the first part of the study and use it as a constant during calculations.
- The geometrical setup of the canal, as depicted in the provided diagram, is essential for accurate calculations.

## License

The code in this repository is provided for educational and research purposes. If you wish to use it in a different context, please ensure proper attribution. 
