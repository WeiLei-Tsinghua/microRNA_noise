# microRNA_noise

Here deposites the code for simulating gene expression noise modulated by microRNA as well as calculating gene expression noise from flow cytometry data.

## Simulate gene expression noise modulated by microRNA

This code is created and tested with MATLAB 2020b.

Run `simulation/solve_noise.m` to simulate noise, and run `simulation/plot_noise.m` to plot the result. An example of the result with default parameters is shown in `simulation/sample_output`.

### Set the type of model

There are three types of preset miRNA regulation models, which can be set by the parameter `type` in Line 7 of `simulation/plot_noise.m`.
For competing RNAs, `type = 1`; for repetitive targets of same miRNAs, `type = 2`; for multiple targets of different miRNAs, `type = 3`. 
 
### Set parameters

All parameters involved in the model should be set before simulation in the `set fixed parameters` section. All parameters that should be altered for simulating different conditions should be set in the `set altered parameters` section.

The simulation step size and range should be assigned in the `set simulation parameters` section. The simulation range determines the range that the gene expression is simulated in by setting the production rate of the target gene’s mRNA (`kT`). We recommend users to set a wide range with a large step size to quickly determine the lower and upper limit of the range, and then narrow the step size to gain a refined simulation result.


## Calculate gene expression noise from flow cytometry data

This code is created and tested with R 3.6.1. ggplot2 is necessary for visualizing the result.

Run `flow_cytometry/plot.r` to calculate noise from the sample data in `flow_cytometry/data`. Sample outputs are shown in `flow_cytometry/result`.

## Contact

For any issues, please contact `weilei92@tsinghua.edu.cn`.
