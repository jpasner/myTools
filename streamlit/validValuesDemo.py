import os
import streamlit as st
import pandas as pd
from snowflake.connector import connect
from dotenv import load_dotenv
from snowflake.connector.pandas_tools import write_pandas
import datetime

# Load credentials from .env file
load_dotenv('.env')

# Connect to Snowflake
conn = connect(
    user=os.getenv('SNOWFLAKE_USER'),
    password=os.getenv('SNOWFLAKE_PASSWORD'),
    account=os.getenv('SNOWFLAKE_ACCOUNT'),
    warehouse=os.getenv('SNOWFLAKE_WAREHOUSE'),
    database="BROTIME",
    schema="LIVINGROOM"
)

# Read valid locations table from Snowflake
locations_query = "SELECT * FROM valid_locations"
locations_df = pd.read_sql(locations_query, conn)

# Read valid folder number table from Snowflake
folder_numbers_query = "SELECT * FROM valid_folder_numbers"
folder_numbers_df = pd.read_sql(folder_numbers_query, conn)

# Display valid locations and folder numbers side by side
st.subheader("Valid Locations and Folder Numbers")
st.write(pd.concat([locations_df, folder_numbers_df], axis=1))
upload_df = pd.DataFrame(columns=['Folder Number', 'Location', 'Timestamp'])

# Streamlit app
st.title("Upload Data")


# Input fields for folder number and location
folder_number = "A0000000001"  # Default value for folder number
folder_number = st.text_input("Folder Number", value=folder_number)  # Text input field for folder number

location = "One"  # Default value for location
location = st.text_input("Location", value=location)  # Text input field for location


# Upload button
if st.button("Upload"):
    # Input fields into new data frame and store current Date and Time in Timestamp column
    new_data = pd.DataFrame({'Folder Number': [folder_number], 'Location': [location], 'Timestamp': [datetime.datetime.now()]})

    # Check if folder number and location are valid
    if folder_number in folder_numbers_df['FOLDER_NUMBER'].values and location in locations_df['LOCATIONS'].values:
        # Concatinate new_data values to upload dataframe
        upload_df = pd.concat([upload_df, new_data], ignore_index=True)

        # Write upload dataframe to Snowflake using write_pandas command
        write_pandas(conn, upload_df, 'Upload_Table', auto_create_table=True)
        st.success("Data uploaded successfully!")
    else:
        st.error("Invalid folder number or location!")

# Display the upload dataframe
st.dataframe(upload_df)
