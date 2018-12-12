import pandas as pd
import numpy as np

xlsx = 'input_data - Copy.xlsx'
output = 'structs-table.rst.txt'

df = pd.read_excel(xlsx, sheet_name='inputdata', skiprows=5)

print(dir(df))
# print(df['Input_for_SCOPE'])

# print(df.iloc[:, 1])
# print(df[df.columns[1]])

with open(output, 'w') as out:
    for _, row in df.iterrows():
        print(row['Variable'])
        if row['Variable'] is np.nan:  # because strings are compared here
            continue
        elif np.isnan(row['Values']):  # because only floats are expected here
            print('here')
            out.write(f"\n{row['Variable']}\n")
            continue
#         out.write(
# f"""
# :{row['Variable']}: {row['Description']}
#
#     :units: {row['Unit']}
#     :type: double
#     :default: {row['Values']}
# """
#         )
        out.write(
    f"""
    * - **{row['Variable']}**
      - {row['Unit']}
      - double
      - {row['Values']}
      - {row['Description']}
    """
        )
