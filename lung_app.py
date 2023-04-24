import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle
import math

# Molecular descriptor calculator
def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')

# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building
def build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('lung_model.pkl','rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    pic50 = prediction.astype('float64')
    
    def cal_ic50(pic50):
        ic50 = (10 **(6-pic50))
        return ic50 
    ic50 = cal_ic50(pic50)

    
    
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    ic50mm_output = pd.Series(ic50, name='IC50 (microM)')
    #activity = pd.Series(['Active' if i < 10 else 'Inactiive' if i > 100 else 'Intermediate' for i in ic50],name="Activity")

    
    df = pd.concat([molecule_name, prediction_output, ic50mm_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)

# Logo image
image = Image.open('lung.png')

st.image(image, use_column_width=True)

# Page title
st.markdown("""
# Bioactivity Prediction App (Chalcones)

This app allows you to predict the bioactivity towards Chalcone compound. `Chalcones` is a drug target for cancer disease.

**Credits**
- App built in `Python` 
- Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Read the Paper]](https://doi.org/10.1002/jcc.21707).
---
""")

# Sidebar
with st.sidebar.title("Upload your data"):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    #uploaded_text = st.text_input("Paste your SMILES")

if st.sidebar.button('Predict'):
    load_data = pd.read_table(uploaded_file, sep=' ', header=None)
    #load_data = read_text(uploaded_text)
    load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)

    st.header('**Original input data**')
    st.write(load_data)

    with st.spinner("Calculating descriptors..."):
        desc_calc()

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors**')
    desc = pd.read_csv('descriptors_output.csv')
    st.write(desc)
    st.write(desc.shape)

    # Read descriptor list used in previously built model
    st.header('**Subset of descriptors from previously built models**')
    Xlist = list(pd.read_csv('descriptor_list.csv').columns)
    desc_subset = desc[Xlist]
    st.write(desc_subset)
    st.write(desc_subset.shape)

    # Apply trained model to make prediction on query compounds
    build_model(desc_subset)
else:
    st.info('Upload input data in the sidebar to start!')
