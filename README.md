In order to build OpenMC and generate example cases, please `cd` into
`fet-stationarity/scripts` and run `python3 setup_project.py -b`. This
will allow you to run OpenMC individually on each example case to generate
output data for data analysis

Since these example cases require significant computational resources to run,
the output data generation step can be skipped altogether by downloading pre-generated
output from a server. In order to build OpenMC, generate example cases, and download
all output data required for running FET stationarity data analysis, please `cd` into 
`fet-stationarity/scripts` and run `python3 setup_project.py -o` instead.

Once all output data has been generated, an example Jupyter notebook that outlines
how all plots found in the paper can be generated. The notebook is located in
`fet-stationarity/examples/fet-plots.ipynb`.
