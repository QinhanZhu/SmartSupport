#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 13:58:32 2024

@author: qinhan
"""

import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from GUI_functions import save_data, load_data, extend_time_impact, find_location
import matplotlib
import pandas as pd



def generate_output():

    total_asset_value = float(Total_asset_value_entry_field.get())    

    Y = []
    for entry in haz_entry_fields[1:8]:
        i = float(entry.get())*total_asset_value
        Y.append(i)
    
    X = [5, 10, 20, 50, 100, 500, 1000]
                 

    # Perform linear interpolation on Y to get the extended Y data series

    extended_X, extended_Y = extend_time_impact(X, Y)
    
    available_financing_resources = 0
    for entry in finance_entry_fields:
        i = float(entry.get())
        available_financing_resources += i
    
    
    financing_gap_year = find_location(extended_Y, available_financing_resources)
    
    # Generate plot
    EFC_ax.clear()
    
    EFC_ax.plot(X,Y)
    
    EFC_ax.set_xlabel('Year events return period (years)')
    
    EFC_ax.set_ylabel('Total damage to your assets (million $)')
    
    vertical_line_financing_gap = financing_gap_year
    
    EFC_ax.axvline(x=vertical_line_financing_gap, color='r', linestyle='--')

    # Add a text box next to the vertical line
    EFC_ax.text(vertical_line_financing_gap + 50, Y[4], f"Financing gap happens every \n{int(financing_gap_year)} years in your country")

    # Add labels and title
    plt.xlabel('X')
    plt.ylabel('Y')
    
    # Update canvas with new plot
    hazard_canvas.draw()
    
    # amount of each financing resource

    values = [None]*5

    values[0] = float(finance_entry_fields[0].get())
    values[1] = float(finance_entry_fields[1].get())
    values[2] = float(finance_entry_fields[2].get())
    values[3] = float(finance_entry_fields[3].get())
    values[4] = float(finance_entry_fields[4].get())

    # find the upper limit of damages each financing resource could cover
    upper_limit = [None]*5

    upper_limit[0] = find_location(extended_Y, values[0])

    upper_limit[1] = find_location(extended_Y, values[0]+values[1])

    upper_limit[2] = find_location(extended_Y, values[0]+values[1]+values[2])

    upper_limit[3] = find_location(extended_Y, values[0]+values[1]+values[2]+values[3])

    upper_limit[4] = find_location(extended_Y, values[0]+values[1]+values[2]+values[3]+values[4])

    # generate the series of usages of each financing resource and the occurance of fiscal gap
    
    assistance = [0]*len(extended_Y)
    diversion = [0] *len(extended_Y)
    taxation = [0] *len(extended_Y)
    domestic = [0] *len(extended_Y)
    MFI = [0] *len(extended_Y)
    gap = [0] *len(extended_Y)


    assistance[:upper_limit[0]] = extended_Y[:upper_limit[0]]
    assistance[upper_limit[0]:] = [values[0]]*len(assistance[upper_limit[0]:])

    diversion[upper_limit[0]:upper_limit[1]] = [x-values[0] for x in extended_Y[upper_limit[0]:upper_limit[1]]]
    diversion[upper_limit[1]:] = [values[1]]*len(diversion[upper_limit[1]:])

    taxation[upper_limit[1]:upper_limit[2]] = [x-values[0]- values[1] for x in extended_Y[upper_limit[1]:upper_limit[2]]]
    taxation[upper_limit[2]:] = [values[2]]*len(taxation[upper_limit[2]:])

    domestic[upper_limit[2]:upper_limit[3]] = [x-values[0]- values[1]-values[2] for x in extended_Y[upper_limit[2]:upper_limit[3]]]
    domestic[upper_limit[3]:] = [values[3]]*len(domestic[upper_limit[3]:])

    MFI[upper_limit[3]:upper_limit[4]] = [x-values[0]- values[1]- values[2]- values[3] for x in extended_Y[upper_limit[3]:upper_limit[4]]]
    MFI[upper_limit[4]:] = [values[4]]*len(MFI[upper_limit[4]:])

    gap[upper_limit[4]:] = [x-values[0]- values[1]- values[2]- values[3] - values[4] for x in extended_Y[upper_limit[4]:]]    

    time_frame = 200 #time span of financing gap plot

    # Generate plot
    
        # gap plot
    gap_ax.clear() # Clear previous plot
    finance_labels = [i.cget('text') for i in finance_entry_label_fields] +['gap']
    col = {finance_entry_label_texts[0]: '#386641', 
          finance_entry_label_texts[1]: '#6A994E', 
          finance_entry_label_texts[2]: '#A7C957',
          finance_entry_label_texts[3]: '#782832',
          finance_entry_label_texts[4]: '#C9182C',
          'gap' : '#FBF7EF'}
    gap_ax.stackplot(extended_X[:time_frame],
                     (assistance[:time_frame],
                      diversion[:time_frame],
                      taxation[:time_frame], 
                      domestic[:time_frame],
                      MFI[:time_frame], 
                      gap[:time_frame]),
                     baseline='zero', 
                     labels = finance_labels,
                     colors = [col.get(l, '#9b59b6') for l in finance_labels])
    gap_ax.set_xlabel('Year event')
    gap_ax.set_ylabel('Values (million $)')
    handles, labels = gap_ax.get_legend_handles_labels()   #get the handles
    gap_ax.legend(reversed(handles), reversed(labels),loc='center left', bbox_to_anchor=(1.05, 0.5),
          fancybox=True, shadow=True, ncol=1, labelspacing = 1.8, fontsize="12")
        # Update canvas with new plot
    gap_canvas.draw()

        #composition plot
    composition_ax.clear()
    indices = [int(20-int(extended_X[0])), int(50-int(extended_X[0])), int(100-int(extended_X[0]))]
    specific_assistance = []
    specific_diversion = []
    specific_taxation = []
    specific_domestic = []
    specific_MFI = []
    specific_gap = []
    
    for i in indices:
        specific_assistance.append(assistance[i])
    arr_specific_assistance = np.array(specific_assistance)
    for i in indices:
        specific_diversion.append(diversion[i])
    arr_specific_diversion = np.array(specific_diversion)
    for i in indices:
        specific_taxation.append(taxation[i])
    arr_specific_taxation = np.array(specific_taxation)
    for i in indices:
        specific_domestic.append(domestic[i])
    arr_specific_domestic = np.array(specific_domestic)
    for i in indices:
        specific_MFI.append(MFI[i])
    arr_specific_MFI = np.array(specific_MFI)
    for i in indices:
        specific_gap.append(gap[i])
    arr_specific_gap = np.array(specific_gap)
    
    
    specific_year_composition = np.concatenate((arr_specific_assistance,
                                                arr_specific_diversion,
                                                arr_specific_taxation,
                                                arr_specific_domestic,
                                                arr_specific_MFI,
                                                arr_specific_gap))
    composition_ax.bar(['20 year','50 year','100 year'],specific_assistance,label=finance_labels[0],color=col.get(finance_labels[0], '#9b59b6'))
    composition_ax.bar(['20 year','50 year','100 year'],specific_diversion,bottom=np.sum(specific_year_composition[:1]),label=finance_labels[1],color=col.get(finance_labels[1], '#9b59b6'))
    composition_ax.bar(['20 year','50 year','100 year'],specific_taxation,bottom=np.sum(specific_year_composition[:2]),label=finance_labels[2],color=col.get(finance_labels[2], '#9b59b6'))
    composition_ax.bar(['20 year','50 year','100 year'],specific_domestic,bottom=np.sum(specific_year_composition[:3]),label=finance_labels[3],color=col.get(finance_labels[3], '#9b59b6'))
    composition_ax.bar(['20 year','50 year','100 year'],specific_MFI,bottom=np.sum(specific_year_composition[:4]),label=finance_labels[4],color=col.get(finance_labels[4], '#9b59b6'))
    composition_ax.bar(['20 year','50 year','100 year'],specific_gap,bottom=np.sum(specific_year_composition[:5]),label=finance_labels[5],color=col.get(finance_labels[5], '#9b59b6'))
    
    handles, labels = composition_ax.get_legend_handles_labels()
    composition_ax.set_xlabel('Selected year events')
    composition_ax.set_ylabel('Value of financing resources or gap (million $)')
    composition_ax.legend(reversed(handles), reversed(labels),loc='center left', bbox_to_anchor=(1.05, 0.5),
          fancybox=True, shadow=True, ncol=1, labelspacing = 1.8,fontsize='12')
    composition_canvas.draw()
    
    #TODO: be aware that selected year should also be larger than X[0]
    
    selected_year = float(selected_year_entry_field.get())
    index = int(selected_year-int(extended_X[0]))
    
    fiscal_gap_event.config(state="normal")  # Enable editing of the output box
    fiscal_gap_event.delete(0, tk.END)  # Clear the current content
    fiscal_gap_event.insert(0, financing_gap_year)  # Insert the new value
    fiscal_gap_event.config(state="readonly")  # Disable editing of the output box
    
    
    selected_assiance = assistance[index]
    selected_diversion = diversion[index]
    selected_taxation = taxation[index]
    selected_domestic = domestic[index]
    selected_MFI = MFI[index]
    selected_gap = gap[index]
      
    finance_resource1.config(state="normal")  # Enable editing of the output box
    finance_resource1.delete(0, tk.END)  # Clear the current content
    finance_resource1.insert(0, selected_assiance)  # Insert the new value
    finance_resource1.config(state="readonly")  # Disable editing of the output box
    
    finance_resource2.config(state="normal")  # Enable editing of the output box
    finance_resource2.delete(0, tk.END)  # Clear the current content
    finance_resource2.insert(0, selected_diversion)  # Insert the new value
    finance_resource2.config(state="readonly")  # Disable editing of the output box
    
    finance_resource3.config(state="normal")  # Enable editing of the output box
    finance_resource3.delete(0, tk.END)  # Clear the current content
    finance_resource3.insert(0, selected_taxation)  # Insert the new value
    finance_resource3.config(state="readonly")  # Disable editing of the output box
     
    finance_resource4.config(state="normal")  # Enable editing of the output box
    finance_resource4.delete(0, tk.END)  # Clear the current content
    finance_resource4.insert(0, selected_domestic)  # Insert the new value
    finance_resource4.config(state="readonly")  # Disable editing of the output box
     
    finance_resource5.config(state="normal")  # Enable editing of the output box
    finance_resource5.delete(0, tk.END)  # Clear the current content
    finance_resource5.insert(0, selected_MFI)  # Insert the new value
    finance_resource5.config(state="readonly")  # Disable editing of the output box
     
    finance_gap.config(state="normal")  # Enable editing of the output box
    finance_gap.delete(0, tk.END)  # Clear the current content
    finance_gap.insert(0, selected_gap)  # Insert the new value
    finance_gap.config(state="readonly")  # Disable editing of the output box
    
    
#TODO: expand the clear range
def clear_entries_plot():
    # Clear entry fields
    for entry in haz_entry_fields:
        entry.delete(0,tk.END)
    
    # Clear plot
    EFC_ax.clear()
    EFC_ax.axis('off')
    hazard_canvas.draw()

    # Clear entry fields
    for entry in finance_entry_fields:
        entry.delete(0,tk.END)
    
    Total_asset_value_entry_field.delete(0,tk.END)
    
    selected_year_entry_field.delete(0,tk.END)
    
    # Clear plot
    gap_ax.clear()
    gap_ax.axis('off')
    gap_canvas.draw()
    
    composition_ax.clear()
    composition_ax.axis('off')
    composition_canvas.draw()

def move_entry_up(entry_index):
    # Move financing resources up
    if entry_index > 0:
        finance_entry_fields[entry_index], finance_entry_fields[entry_index - 1] = finance_entry_fields[entry_index - 1], finance_entry_fields[entry_index]
        finance_entry_label_fields[entry_index], finance_entry_label_fields[entry_index - 1] = finance_entry_label_fields[entry_index - 1], \
        finance_entry_label_fields[entry_index]
        refresh_entry_labels()

def move_entry_down(entry_index):
    # Move financing resources down
    if entry_index < len(finance_entry_fields) - 1:
        finance_entry_fields[entry_index], finance_entry_fields[entry_index + 1] = finance_entry_fields[entry_index + 1], finance_entry_fields[entry_index]
        finance_entry_label_fields[entry_index], finance_entry_label_fields[entry_index + 1] = finance_entry_label_fields[entry_index + 1], \
        finance_entry_label_fields[entry_index]
        refresh_entry_labels()

def refresh_entry_labels():
    # Refresh the labels for the financing resource
    for i, entry in enumerate(finance_entry_fields):
        entry.grid(row=i+10, column=1, padx=10, pady=5)
        if i > 0:
            move_up_buttons[i - 1].grid(row=i+10, column=2, padx=5, pady=2)
        if i < len(finance_entry_fields) - 1:
            move_down_buttons[i].grid(row=i+10, column=3, padx=5, pady=2)

    for i, entry_label in enumerate(finance_entry_label_fields):
        entry_label.grid(row=i+10, column=0)




# Create the main window and define design parameters

window = tk.Tk()
window.title("Final GUI")

# Apply theme
window.tk.call("source", "azure.tcl")
window.tk.call("set_theme", "dark")

# Font of the plots
font = {'family' : 'Georgia',
        'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)


# Create frame for input entries

Input_FrameLabel = ttk.Label(text="Define your risk profile and \navailable financing resources", font=("Georgia", 24))

Input_frame = ttk.LabelFrame(window,labelwidget=Input_FrameLabel, padding=(20, 10))
Input_frame.grid(row=0, column=0, padx=10, pady=10,sticky='nsew')

hazard_input_title_label = ttk.Label(Input_frame, text="Define the risk profile of your country", 
                                     font=("Georgia", 18, "bold"))
hazard_input_title_label.grid(row=0, column=0, columnspan=2, padx = 10, pady = 10)

    # Define the entries for hazard
haz_entry_fields = []
haz_entry_label_texts = ['Probability of first damage','% of asset damage by 5 year event','% of asset damage by 10 year event','% of asset damage by 20 year event','% of asset damage by 50 year event','% of asset damage by 100 year event','% of asset damage by 500 year event','% of asset damage by 1000 year event']
haz_entry_label_fields = []

for i in range(8):
    haz_entry = ttk.Entry(Input_frame)
    haz_entry.grid(row=i+1, column=1, padx=10, pady=5, sticky="w")
    haz_entry_fields.append(haz_entry)
    entry_label = ttk.Label(Input_frame,text = haz_entry_label_texts[i],font=("Georgia", 18))
    entry_label.grid(row=i+1, column=0, padx=10, pady=5, sticky="w")
    haz_entry_label_fields.append(entry_label)



finance_title_label = ttk.Label(Input_frame, text="Define the available financing resources (unit: million $)", 
                                font=("Georgia", 18, "bold"))
finance_title_label.grid(row=9, column=0, padx=10, pady=10, columnspan=2)

    # Create entries for available financing resources
    
finance_entry_fields = []
finance_entry_label_texts = ['Assistance','Diversion','Taxation','Domestic\ncredit','MFI credit \nBond credit']
finance_entry_label_fields = []
move_up_buttons = []
move_down_buttons = []
for i in range(5):
    entry = ttk.Entry(Input_frame)
    entry.grid(row=i+10, column=1, padx=10, pady=5)
    finance_entry_fields.append(entry)
    entry_label = ttk.Label(Input_frame,text = finance_entry_label_texts[i], font=("Georgia", 18))
    entry_label.grid(row=i+10, column=0, padx=10, pady=5, sticky="w")
    finance_entry_label_fields.append(entry_label)

    if i > 0:
        move_up_button = ttk.Button(Input_frame, text="↑", width=1, command=lambda index=i: move_entry_up(index))
        move_up_buttons.append(move_up_button)
        move_up_button.grid(row=i+10, column=2, padx=5, pady=2)

    if i < 4:
        move_down_button = ttk.Button(Input_frame, text="↓", width=1, command=lambda index=i: move_entry_down(index))
        move_down_buttons.append(move_down_button)
        move_down_button.grid(row=i+10, column=3, padx=5, pady=2)

Total_asset_value_entry_field = ttk.Entry(Input_frame)
Total_asset_value_entry_field.grid(row=15, column=1, padx=10, pady=5)
Total_asset_value_entry_label_field = ttk.Label(Input_frame,text = 'Value of assets exposed to climate hazard',
                                            font=("Georgia", 18))
Total_asset_value_entry_label_field.grid(row=15, column=0, padx=10, pady=5,sticky="w")

selected_year_entry_field = ttk.Entry(Input_frame)
selected_year_entry_field.grid(row=16, column=1, padx=10, pady=5)
selected_year_entry_label_field = ttk.Label(Input_frame,text = 'Select a specific year event to \nfind out the optimal financing strategy',
                                            font=("Georgia", 18))
selected_year_entry_label_field.grid(row=16, column=0, padx=10, pady=5,sticky="w")


data_fields = haz_entry_fields
for x in finance_entry_fields:
    data_fields.append(x) 
data_fields.append(Total_asset_value_entry_field)
data_fields.append(selected_year_entry_field)


    # Generate button
generate_button = ttk.Button(Input_frame, text="Generate", command=generate_output)
generate_button.grid(row=17, column = 0, padx=10, pady=5, sticky="w")

    # Clear button
clear_button = ttk.Button(Input_frame, text="Clear", command=clear_entries_plot)
clear_button.grid(row=17, column = 1, padx=10,  pady=5, sticky="w")

haz_status_label = ttk.Label(Input_frame)
haz_status_label.grid(row=19, column=0, columnspan=2, padx=10, pady=5)

    # Save button
save_button = ttk.Button(Input_frame, text="Save Data", command=lambda: save_data(data_fields, haz_status_label))
save_button.grid(row=18, column=0, padx=10, pady=5, sticky="w")


    # Load button
load_button = ttk.Button(Input_frame, text="Load Data", command=lambda: load_data(data_fields, haz_status_label))
load_button.grid(row=18, column=1, padx=10, pady=5, sticky="w")



# Create frame for plot

Output_FrameLabel=ttk.Label(text="Results of your risk profile and financing gap", font=("Georgia", 24))
Output_frame = ttk.LabelFrame(window,labelwidget=Output_FrameLabel, padding=(20, 10))
Output_frame.grid(row=0, column=1, padx=10, pady=10, sticky='nsew')



EFC_frame = ttk.Frame(Output_frame, padding=(20, 10))
EFC_frame.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

EFC_frame_label = ttk.Label(EFC_frame, text= 'Risk distribution and occurance of financing gap', font=("Georgia", 20, 'bold'), padding=(20, 10))
EFC_frame_label.pack()

    # Create Matplotlib figure and canvas

EFC_fig, EFC_ax = plt.subplots(figsize=(6, 4))
EFC_ax.axis('off')
hazard_canvas = FigureCanvasTkAgg(EFC_fig, master=EFC_frame)
hazard_canvas.draw()
EFC_fig.patch.set_alpha(0)
hazard_canvas.get_tk_widget().pack()


# Create frame for financing gap plot
gap_frame = ttk.Frame(Output_frame)
gap_frame.grid(row=1, column=0, padx=10, pady=10, sticky='nsew')

gap_frame_label = ttk.Label(gap_frame, text= 'Amount of resources deployed to \ncover the financing needs', font=("Georgia", 20, 'bold'), padding=(20, 10))
gap_frame_label.pack()


    # Create Matplotlib figure and canvas

gap_fig, gap_ax = plt.subplots(figsize=(6, 4))
gap_ax.axis('off')
box = gap_ax.get_position()
gap_ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
gap_canvas = FigureCanvasTkAgg(gap_fig, master=gap_frame)
gap_canvas.draw()
gap_canvas.get_tk_widget().pack()


# Create frame for financing composition plot for specific year events
composition_frame = ttk.Frame(Output_frame)
composition_frame.grid(row=1, column=1, padx=10, pady=10,sticky='nsew')

composition_frame_label = ttk.Label(composition_frame, text= 'Finance resources deployed to \ncover the 20, 50, and 100 year event', font=("Georgia", 20, 'bold'), padding=(20, 10))
composition_frame_label.pack()

    # Create Matplotlib figure and canvas

composition_fig, composition_ax = plt.subplots(figsize=(6, 4))
composition_ax.axis('off')
box = composition_ax.get_position()
composition_ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
composition_canvas = FigureCanvasTkAgg(composition_fig, master=composition_frame)
composition_canvas.draw()
composition_canvas.get_tk_widget().pack()

#TODO: make a loop
# Create frame for quantative data report

data_frame = ttk.Frame(Output_frame)
data_frame.grid(row=0, column = 1, padx = 10, pady = 10,sticky='nsew')

data_frame_label = ttk.Label(data_frame, text= 'Finance resources deployed to \ncover your specified year event', font=("Georgia", 20, 'bold'), padding=(20, 10))
data_frame_label.grid(row=0, column=0, columnspan=2, padx = 10, pady =5)

fiscal_gap_event_label = ttk.Label(data_frame, text='Year event triggering fiscal gap:', font=("Georgia", 18))
fiscal_gap_event_label.grid(row=1, column=0, padx=10, pady=5, sticky = 'w')

fiscal_gap_event = ttk.Entry(data_frame, state='readonly')
fiscal_gap_event.grid(row=1, column=1, padx=10, pady=5)

selected_year = selected_year_entry_field.get()

composition_label = ttk.Label(data_frame, text='To cover the financing need for the intended year event, \nfrom the following resources you have', font=("Georgia", 18))
composition_label.grid(row=2, padx=10, columnspan=2, pady=5, sticky='w')


finance_labels = [i.cget('text') for i in finance_entry_label_fields] +['gap']


finance_resource1_label = ttk.Label(data_frame, text=finance_labels[0], font=("Georgia", 18))
finance_resource1_label.grid(row=3, column=0, padx=10, pady=5, sticky = 'w')

finance_resource1 = ttk.Entry(data_frame, state='readonly')
finance_resource1.grid(row=3, column=1, padx=10, pady=5)

finance_resource2_label = ttk.Label(data_frame, text=finance_labels[1], font=("Georgia", 18))
finance_resource2_label.grid(row=4, column=0, padx=10, pady=5, sticky = 'w')

finance_resource2 = ttk.Entry(data_frame, state='readonly')
finance_resource2.grid(row=4, column=1, padx=10, pady=5)

finance_resource3_label = ttk.Label(data_frame, text=finance_labels[2], font=("Georgia", 18))
finance_resource3_label.grid(row=5, column=0, padx=10, pady=5, sticky = 'w')

finance_resource3 = ttk.Entry(data_frame, state='readonly')
finance_resource3.grid(row=5, column=1, padx=10, pady=5)

finance_resource4_label = ttk.Label(data_frame, text=finance_labels[3], font=("Georgia", 18))
finance_resource4_label.grid(row=6, column=0, padx=10, pady=5, sticky = 'w')

finance_resource4 = ttk.Entry(data_frame, state='readonly')
finance_resource4.grid(row=6, column=1, padx=10, pady=5)

finance_resource5_label = ttk.Label(data_frame, text=finance_labels[4], font=("Georgia", 18))
finance_resource5_label.grid(row=7, column=0, padx=10, pady=5, sticky = 'w')

finance_resource5 = ttk.Entry(data_frame, state='readonly')
finance_resource5.grid(row=7, column=1, padx=10, pady=5)

gap_label = ttk.Label(data_frame, text='The finance gap would be', font=("Georgia", 18))
gap_label.grid(row=8, column=0, padx=10, pady=5, sticky = 'w')

finance_gap = ttk.Entry(data_frame, state='readonly')
finance_gap.grid(row=8, column=1, padx=10, pady=5)


window.update()
# window.minsize(window.winfo_width(), window.winfo_height())
window.maxsize(2000, 2000)
x_cordinate = int((window.winfo_screenwidth() / 2) - (window.winfo_width() / 2))
y_cordinate = int((window.winfo_screenheight() / 2) - (window.winfo_height() / 2))
window.geometry("+{}+{}".format(x_cordinate, y_cordinate-20))
# Run the GUI
window.mainloop()