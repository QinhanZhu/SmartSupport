import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import bisect

def save_data(entries, status_label):
    data = [entry.get() for entry in entries]
    file_path = filedialog.asksaveasfilename(filetypes=(("Text Files", "*.txt"), ("All Files", "*.*")))
    if file_path:
        with open(file_path, "w") as file:
            for item in data:
                file.write(item + "\n")
        status_label.config(text="Data saved successfully!")

def load_data(entries, status_label):
    file_path = filedialog.askopenfilename(filetypes=(("Text Files", "*.txt"), ("All Files", "*.*")))
    if file_path:
        data = []
        try:
            with open(file_path, "r") as file:
                for line in file:
                    data.append(line.strip())
            if len(data) >= 0:
                for i in range(len(entries)):
                    entries[i].delete(0, tk.END)
                    entries[i].insert(tk.END, data[i])
                status_label.config(text="Data loaded successfully!")
            else:
                status_label.config(text="Invalid data format!")
        except FileNotFoundError:
            status_label.config(text="File not found!")


def extend_time_impact(X,Y):
    
    extended_Y = []
    extended_X = []

    if X[0] == int(X[0]):
        for i in range(len(Y) - 1):
            y_current, y_next = Y[i], Y[i + 1]
            slope = (y_next - y_current) / (X[i + 1] - X[i])
            interpolated_points = [y_current + slope * (x - X[i]) for x in range(int(X[i]), X[i + 1])]
            extended_Y.extend(interpolated_points)
        extended_Y.append(Y[-1])  # Append the last point
        for i in range(len(X) - 1):
            x_current, x_next = X[i], X[i + 1]
            dx = x_next - x_current
            num_points = int(dx)
            interpolated_points = np.linspace(x_current, x_next, num_points + 1)
            extended_X.extend(interpolated_points[:-1])
        extended_X.append(X[-1])
        
    elif int(X[0])+1 == X[1]:
        for i in range(1,len(Y) - 1):
            y_current, y_next = Y[i], Y[i + 1]
            slope = (y_next - y_current) / (X[i + 1] - X[i])
            interpolated_points = [y_current + slope * (x - X[i]) for x in range(X[i], X[i + 1])]
            extended_Y.extend(interpolated_points)
        extended_Y.append(Y[-1])  # Append the last point
        extended_Y.insert(0, Y[0])
        for i in range(len(X) - 1):
            x_current, x_next = X[i], X[i + 1]
            dx = x_next - x_current
            num_points = int(dx)
            interpolated_points = np.linspace(x_current, x_next, num_points + 1)
            extended_X.extend(interpolated_points[:-1])
        extended_X.append(X[-1])
        extended_X.insert(0, X[0])
        
    else:
        for i in range(1,len(Y) - 1):
            y_current, y_next = Y[i], Y[i + 1]
            slope = (y_next - y_current) / (X[i + 1] - X[i])
            interpolated_points = [y_current + slope * (x - X[i]) for x in range(X[i], X[i + 1])]
            extended_Y.extend(interpolated_points)
        extended_Y.append(Y[-1])  # Append the last point
        first_slope = (Y[1] - Y[0])/ (X[1] - X[0])
        interpolated_points = [Y[0] + first_slope * (x - X[0]) for x in range(int(X[0])+1, X[1])]
        extended_Y = interpolated_points + extended_Y
        extended_Y.insert(0, Y[0])
        
        for i in range(1,len(X) - 1):
            x_current, x_next = X[i], X[i + 1]
            dx = x_next - x_current
            num_points = int(dx)
            interpolated_points = np.linspace(x_current, x_next, num_points + 1)
            extended_X.extend(interpolated_points[:-1])
        extended_X.append(X[-1])
        len_first_xs = int(X[1] - X[0])
        first_xs = [X[0]]
        for i in range(len_first_xs):
            first_xs.append(int(X[0])+1+i)
        extended_X = first_xs + extended_X
        
    return extended_X, extended_Y

def find_location(series,threshold):
    
    left, right = 0, len(series) - 1
    while left <= right:
        mid = left + (right - left) // 2
        if series[mid] == threshold:
            return mid+1  # Value x found at index mid
        elif series[mid] < threshold:
            left = mid + 1
        else:
            right = mid - 1
    # # If value x is not found, return the location where x would be inserted
    return right

def calc_fiscal_gap(year, X, assistance, diversion, taxation, dom_credit, MFI_credit):
    
    damage = X[year]
    max_fr_assitance = assistance * damage
    max_fr_diversion = diversion
    max_taxation = taxation
    max_domcre = dom_credit
    max_MFIcre = MFI_credit

    