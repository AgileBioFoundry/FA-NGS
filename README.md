# Fluorescent Amplification of Next Generation Sequencing Libraries Tools (FA-NGS Tools)

This is a software utility for calculating pooling volumes and graphing quality control data of next generation sequencing (NGS) libraries amplified with fluorescent intercalating agents and monitored by real-time quantitative PCR (qPCR). The tool is used to calculate pooling volumes of individual NGS libraries based on experimentally determined fluorescence values. To do this, fluorescent values are used as a proxy for relative concentration to determine the ratio of transfer volumes. The output is a “transfer file” that specifies the transfer volume  and destination well of a given library. The tool can also be used to assess the quality of individual next generation sequencing libraries by melting curve analysis. The output are small-multiples graphs of the melting curves of individual sequencing libraries according to the specified plate layout.

## Getting Started
Use the example files and Jupyter notebooks to quickly learn how to use the tools

### Prerequisites

* Numpy
* Pandas
* Matplotlib
* Seaborn
* Streamlit
* Base64

# Authors
Megan E. Garber

# Acknowledgments 
This project was funded by the Ecosystems and Networks Integrated with Genes and Molecular Assemblies a Scientific Focus Area Program at Lawrence Berkeley National Laboratory and is supported by the U.S. Department of Energy, Office of Science, Office of Biological & Environmental Research under contract number DE-AC02–05CH11231 between Lawrence Berkeley National Laboratory and the U. S. Department of Energy.

# Licencing
FA-NGS tools Copyright (c) 2019, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

 

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov referring to " FA-NGS tools" (LBNL Ref 2019-158)."

 

NOTICE.  This software was developed under funding from the U.S. Department of Energy.  As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, and perform publicly and display publicly.  The U.S. Government is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

