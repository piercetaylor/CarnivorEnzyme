"""CarnivorEnzyme Streamlit app entry point."""
import streamlit as st

st.set_page_config(page_title="CarnivorEnzyme Atlas", layout="wide")
st.title("CarnivorEnzyme: Structural Atlas")
st.info("Run the Snakemake pipeline first to populate the database.")
