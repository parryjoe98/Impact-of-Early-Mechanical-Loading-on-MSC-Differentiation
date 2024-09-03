import pandas as pd
import matplotlib.pyplot as plt

Load = ['1e5', '25e4','5e5', '1e6']

iteration = []

for loads in Load:
    
    Loading = loads

    for day in range(1,8):
        Loading_day = day
    
        
        # Define the base file path and the number of files
        base_path = f'C:/Users/joepa/OneDrive - Imperial College London/Desktop/Project full backups/Full/Day{Loading_day}/Load{Loading}/'
        file_prefix = 'E_elem_'
        file_suffix = '.xls'
        num_files = 35
        
        # Define the threshold value to look for in column Q
        threshold_value = 2e9  # replace with your threshold value
        
        first_found = False
        
        # Function to add specific rows to Elem_row
        def add_specific_rows(df, base_row, additional_row):
            row_base = df.iloc[base_row - 1:base_row]  # Subtract 1 for zero-based indexing
            row_additional = df.iloc[additional_row - 1:additional_row]
            return pd.concat([row_base, row_additional], ignore_index=True)
        
        for i in range(1, num_files + 1):
            # Generate the file path for the current file
            file_path = f'{base_path}{file_prefix}{i}{file_suffix}'
            
            # Load the Excel file
            df = pd.read_excel(file_path)
            
            # Extract specific rows by index (each row 8 long)
            Elem_row_1 = df.iloc[3044:3052]
            Elem_row_2 = df.iloc[3052:3060]
            Elem_row_3 = df.iloc[3060:3068]
            Elem_row_4 = df.iloc[3068:3076]
            Elem_row_5 = df.iloc[3076:3084]
            Elem_row_6 = df.iloc[3084:3092]
            Elem_row_7 = df.iloc[3092:3100]
            Elem_row_8 = df.iloc[3100:3108]
            Elem_row_9 = df.iloc[3108:3116]
            Elem_row_10 = df.iloc[3116:3124]
            Elem_row_11 = df.iloc[3124:3132]
            Elem_row_12 = df.iloc[3132:3140]
            Elem_row_13 = df.iloc[3140:3148]
            
            
            # Add specific rows to Elem_row_1 to Elem_row_13
            Elem_row_1 = pd.concat([Elem_row_1, add_specific_rows(df, 136, 197)], ignore_index=True)
            Elem_row_2 = pd.concat([Elem_row_2, add_specific_rows(df, 137, 198)], ignore_index=True)
            Elem_row_3 = pd.concat([Elem_row_3, add_specific_rows(df, 138, 199)], ignore_index=True)
            Elem_row_4 = pd.concat([Elem_row_4, add_specific_rows(df, 139, 200)], ignore_index=True)
            Elem_row_5 = pd.concat([Elem_row_5, add_specific_rows(df, 140, 201)], ignore_index=True)
            Elem_row_6 = pd.concat([Elem_row_6, add_specific_rows(df, 141, 202)], ignore_index=True)
            Elem_row_7 = pd.concat([Elem_row_7, add_specific_rows(df, 142, 203)], ignore_index=True)
            Elem_row_8 = pd.concat([Elem_row_8, add_specific_rows(df, 143, 204)], ignore_index=True)
            Elem_row_9 = pd.concat([Elem_row_9, add_specific_rows(df, 144, 205)], ignore_index=True)
            Elem_row_10 = pd.concat([Elem_row_10, add_specific_rows(df, 145, 206)], ignore_index=True)
            Elem_row_11 = pd.concat([Elem_row_11, add_specific_rows(df, 146, 207)], ignore_index=True)
            Elem_row_12 = pd.concat([Elem_row_12, add_specific_rows(df, 147, 208)], ignore_index=True)
            Elem_row_13 = pd.concat([Elem_row_13, add_specific_rows(df, 148, 209)], ignore_index=True)
        
            
            # Define the output file path
            output_path = f'C:/Users/joepa/OneDrive - Imperial College London/Desktop/Project full backups/Results/BB/Load{Loading}/extracted_rows_{i}.xlsx'
            
            # Save extracted rows to a new Excel file with multiple sheets
            with pd.ExcelWriter(output_path) as writer:
                Elem_row_1.to_excel(writer, sheet_name='Elem_row_1', index=False)
                Elem_row_2.to_excel(writer, sheet_name='Elem_row_2', index=False)
                Elem_row_3.to_excel(writer, sheet_name='Elem_row_3', index=False)
                Elem_row_4.to_excel(writer, sheet_name='Elem_row_4', index=False)
                Elem_row_5.to_excel(writer, sheet_name='Elem_row_5', index=False)
                Elem_row_6.to_excel(writer, sheet_name='Elem_row_6', index=False)
                Elem_row_7.to_excel(writer, sheet_name='Elem_row_7', index=False)
                Elem_row_8.to_excel(writer, sheet_name='Elem_row_8', index=False)
                Elem_row_9.to_excel(writer, sheet_name='Elem_row_9', index=False)
                Elem_row_10.to_excel(writer, sheet_name='Elem_row_10', index=False)
                Elem_row_11.to_excel(writer, sheet_name='Elem_row_11', index=False)
                Elem_row_12.to_excel(writer, sheet_name='Elem_row_12', index=False)
                Elem_row_13.to_excel(writer, sheet_name='Elem_row_13', index=False)
            
            # Check the extracted rows for the value
            if not first_found:
                # Load the newly created Excel file to check for the value
                xls = pd.ExcelFile(output_path)
                
                # Check each sheet within the file
                for sheet_name in xls.sheet_names:
                    df_extracted = pd.read_excel(output_path, sheet_name=sheet_name)
                    
                    # Check if the dataframe has at least 16 columns
                    if df_extracted.shape[1] >= 17:
                        # Check if all values in the 16th column meet or exceed the threshold
                        if (df_extracted.iloc[:, 17] >= threshold_value).all():
                            print(f'Loading of {Loading}Pa at Day{Loading_day} in iteration {i}')
                            #print(df_extracted.iloc[:,17])
                            #Day_iter = i(
                            iteration.append(i)
                            first_found = True
                            break

plt.figure(figsize=(10, 6))
plt.plot(range(1, 8), iteration, linestyle='-', color='b')
plt.title(f'Day that Bony Bridging occured with Load {Load}')
plt.xlabel('Day of Loading')
plt.ylabel('Day that Bony Bridging occured')
plt.ylim(ymin=0,ymax=30)
plt.show()