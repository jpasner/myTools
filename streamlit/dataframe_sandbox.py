import streamlit as st
import pandas as pd
import datetime
from dotenv import load_dotenv
import os
import snowflake.connector
from snowflake.connector.pandas_tools import write_pandas

# Load Snowflake credentials from .env file
load_dotenv('.env')

# Create an empty DataFrame
df = pd.DataFrame(columns=['Number 1', 'Number 2', 'Timestamp'])

# Streamlit app
st.title("Number Input App")

# Input fields for two numbers
number1 = st.number_input("Enter Number 1", value=0)
number2 = st.number_input("Enter Number 2", value=0)

# Button to store the numbers in the DataFrame and Snowflake database
if st.button("Store Numbers"):
    # Get the current system time
    timestamp = datetime.datetime.now()

    # Create a new DataFrame with the numbers and timestamp
    new_data = pd.DataFrame({'Number 1': [number1], 'Number 2': [number2], 'Timestamp': [timestamp]})

    # Concatenate the new DataFrame with the existing DataFrame
    df = pd.concat([df, new_data], ignore_index=True)

    # Show the DataFrame
    st.dataframe(df)

    # Connect to Snowflake using credentials from .env file
    conn = snowflake.connector.connect(
        user=os.getenv('SNOWFLAKE_USER'),
        password=os.getenv('SNOWFLAKE_PASSWORD'),
        account=os.getenv('SNOWFLAKE_ACCOUNT'),
        warehouse=os.getenv('SNOWFLAKE_WAREHOUSE'),
        schema=os.getenv('SNOWFLAKE_SCHEMA'),
        database=os.getenv('SNOWFLAKE_DATABASE'),
    )

    # Manually set SCHEMA since write_pandas apparently doesn't know how to interpret the snowflake connection object properly
    cursor = conn.cursor()
    cursor.execute("USE SCHEMA FIELD_ONE")

    # Load the DataFrame into Snowflake using the to_sql command
    write_pandas(conn, df, 'numbers', auto_create_table=True)

    # Show success message
    st.success("Successfully stored the numbers and timestamp in the Snowflake database.")

    # Close the cursor and connection
    conn.close()