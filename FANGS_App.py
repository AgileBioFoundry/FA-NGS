"""FA-NGS Tools is user interface built for calculating
    pooling volumes of next generation sequencing libraries amplified with fluorescent dyes"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import base64
from FANGS import *
st.set_page_config(layout="wide")
st.markdown("# Welcome to FA-NGS Tools")

st.markdown("FA-NGS Tools is used to calculate pooling volumes of next generation sequencing libraries amplified with fluorescent dyes.")


#initialize
quad, letters, numbers = None, None, None

# Input plate setup information
st.sidebar.markdown("## Plate Specification")
plate_size = st.sidebar.selectbox("Select a plate option", 
             options = ["384 well","96 well","custom"])

if plate_size == "custom":
    st.sidebar.markdown("### Specify the wells (columns and rows) of your plate:")
    letters = st.sidebar.text_input("Input active plate rows (e.g. ABCDEFG)", "ABCDEFGHIJK")
    numbers = st.sidebar.slider("Drag the slider to indicate the maximum column number",
                        1,24,24,1 )
    
    
    
elif plate_size == "384 well":
    plate_size = 384
    quad_check = st.sidebar.checkbox("Do you want to specify a quadrant?")
    
    if quad_check:
        quad = st.sidebar.slider("Which quadrant is active?",1,4,1,1 )
    else:
        quad = None
else:
    plate_size = 96
    pass


st.markdown("## File Uploads")

col1, col2 = st.beta_columns(2)

#col1 will be END RFU file
col1.markdown("#### Upload the END RFU file")
EndFile = col1.file_uploader("End RFU File")

# col2 will be Melt curve file
col2.markdown("#### Upload the Melt Curve file")
MeltFile = col2.file_uploader("Melt Curve File")
MeltUploadButton = col2.button("Upload Melt File")

@st.cache()
def MeltAnanalysis():
    melt_analysis = Analysis(MeltFile, plate_size= plate_size, quad= quad, 
                 letters= letters, numbers = numbers, failed_wells= [])
    return melt_analysis

@st.cache()
def EndAnalysis():
    end_analysis = Analysis(EndFile, plate_size= plate_size, 
                            quad= quad, letters= letters, numbers = numbers, 
                            failed_wells=[])
    return end_analysis

def PoolSamples():
    Pooling = Pooling_Calculator(EndFile, plate_size= plate_size,
                                     quad= quad, letters= letters, 
                                     numbers = numbers, 
                                     failed_wells=failed_wells)
    return Pooling

def PlotMeltCurve(analysis):
    plot = melt_analysis.plot_melt_curve()
    return plot


def PlotEndRFU(analysis):
    plot = end_analysis.plot_end_RFU()
    return plot


# End RFU and Melt Figures
if MeltFile and MeltUploadButton:
    melt_analysis = MeltAnanalysis()
    fig = PlotMeltCurve(melt_analysis)
    
    col2.write(fig)
else:
    pass


if EndFile:
    end_analysis = EndAnalysis()
    
    col1.write(PlotEndRFU(end_analysis))
    
    
    
#Well removal
well_removal = st.beta_expander("Do you want to remove any wells from the calculation?")
failed_wells = list()
with well_removal:
    number_of_wells = st.number_input("How many wells failed QC?", 
                    min_value =0, max_value =None,  value = 1, step = 1)
    if number_of_wells > 0:
        st.markdown("Input which wells failed:")
        for i in range(number_of_wells):
            wellLabel = st.text_input("", key = i)
            if wellLabel:
                failed_wells.append(wellLabel)

        st.markdown("You are removing wells: {}".format(", ".join(failed_wells)))
    else:
        st.markdown("You are not removing any wells {} from the calculation"\
                    .format(",".join(failed_wells)))



#Echo calculations

PoolingButton = st.beta_expander("Show Pooling Volumes")

Pooling = PoolSamples()

with PoolingButton: 
    if EndFile:
        Pooling_file =Pooling.make_pooling_file()
        
        st.markdown("Echo Pooling File:")
        st.write(Pooling_file)
        
        PoolingPlot = st.checkbox("Show Pooling Plot")
        if PoolingPlot:
            st.write(Pooling.plot_transfer_volumes())
        else:
            pass

    else:
        st.error("No End RFU file uploaded")

#File Download
DownloadButton = st.beta_expander("Download Pooling File")

with DownloadButton:
    if EndFile:

        csv = Pooling_file.to_csv(index=False)
        b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
        href = f'<a href="data:file/csv;base64,{b64}">Download CSV File</a> (right-click and save as &lt;some_name&gt;.csv)'
        st.markdown(href, unsafe_allow_html=True)

    else:
        st.error("No End RFU file uploaded")
