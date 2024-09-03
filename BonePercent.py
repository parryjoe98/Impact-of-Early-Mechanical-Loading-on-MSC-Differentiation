import pandas as pd
import matplotlib.pyplot as plt

Load = ['1e5', '25e4', '5e5', '1e6']

# Initialize a dictionary to store averages for each load and day
averages_dict = {load: [] for load in Load}

day_averages_dict = {day: [] for day in range(1, 8)}

for loads in Load:
    
    Loading = loads

    for day in range(1,8):
        Loading_day = day
    
        
        # Define the base file path and the number of files
        base_path = f'C:/Users/joepa/OneDrive - Imperial College London/Desktop/Project full backups/Full/Day{Loading_day}/Load{Loading}/'
        file_prefix = 'E_elem_'
        file_suffix = '.xls'
        num_files = 35
    
        num_elem = 3392
        
        # Define the threshold value to look for in column Q
        threshold_value = 2e9  # replace with your threshold value
        
        # Initialize a list to store the percentage of immature bone for each day
        percentages = []
        
        growth_rate = []
        
        for i in range(1, num_files + 1):
            # Generate the file path for the current file
            file_path = f'{base_path}{file_prefix}{i}{file_suffix}'
            
            # Load the Excel file
            df = pd.read_excel(file_path)
            
            # Select the 16th column (column Q)
            cell_check = df.iloc[:, 17]
            
            # Initialize the total count of immature bone cells for the current file
            total = 0
            
            # Iterate through each cell value in the 16th column
            for value in cell_check:
                if value > threshold_value:
                    total += 1
            
            ##Calculate the percentage of immature bone
            Percentage_Bone = (total / num_elem) * 100
            
            #print(f'Loading of {Load} on Day{Loading_day} Day {i}: {Percentage_Bone}%')
    
            
            ##Print the percentage for the current day
            
            #if i in (7, 14, 21, 28, 35):
                #print(f'Day of loading {Loading_day} with {Loading} Day {i}: {Percentage_Bone:.3f}%')
                
            #Store the percentage in the list
            percentages.append(Percentage_Bone)
            
            rate_of_change = [(percentages[i] - percentages[i-1]) / percentages[i-1] for i in range(9, len(percentages))] 
            if i == 35:
                #print(f"Day {day} of loading {rate_of_change:.2f}%")
                
                #Average is the average percentage change per iteration/day
                average = sum(rate_of_change) / len(rate_of_change)
                average = average*100
                
                # Append the average to the list of averages for the current load and day
                averages_dict[Loading].append(average)
                
                #print(F'average rate of growth for Day {day} loading is {average:.2f}%')
                #print(f'Loading on Day {day} with {Loading} resulted in a growth rate of {average:.2f}')
                
                # Append the average to the list for the current day across all loads
                day_averages_dict[Loading_day].append(average)
                
                print(f'Loading on Day {Loading_day} with {Loading} resulted in a growth rate of {average:.2f}%')
                
# Calculate the overall average of the averages for each load
overall_averages = {}
for load, averages in averages_dict.items():
    if averages:
        overall_average = sum(averages) / len(averages)
        overall_averages[load] = overall_average
        print(f"Overall average growth rate for load {load}: {overall_average:.2f}%")
        
# Calculate the overall average growth rate for each day across all loads
overall_day_averages = {}
for day, averages in day_averages_dict.items():
    if averages:
        overall_average = sum(averages) / len(averages)
        overall_day_averages[day] = overall_average
        print(f"Overall average growth rate for Day {day} across all loads: {overall_average:.2f}%")


# Print the overall averages dictionary
#print("Overall averages for each load:", overall_averages)
        
                
                            
        # Print the list of percentages for all days
        #print('Percentages of immature bone for each day:', percentages)
        
        # Plot the percentages as a line graph
        plt.figure(figsize=(10, 6))
        plt.plot(range(1, num_files + 1), percentages, linestyle='-', color='b')
        plt.title(f'Percentage of Immature Bone Over 39 Days with Loading of {Load} on Day {Loading_day}')
        plt.xlabel('Day')
        plt.ylabel('Percentage of Immature Bone (%)')
        plt.xticks(range(1, num_files + 1))
        plt.show()