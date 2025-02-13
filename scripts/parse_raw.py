import re
import numpy as np


def parse_amortized_cost(filename):
    # Initialize arrays
    gc_size_array = []
    stack_cost_array = []
    tsc_stack_cost_array = []

    # Read the file
    with open(filename, "r") as file:
        lines = file.readlines()

    # Initialize variables
    T = 0

    # Process each line
    for line in lines:
        if line.startswith("T:"):
            T = int(line.split(":")[1].strip())
        elif line.startswith("gc size:"):
            gc_size = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
            gc_size_array.append(gc_size / T * 1024)
        elif line.startswith("stack cost:"):
            stack_cost = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
            stack_cost_array.append(stack_cost / T * 1024)
        elif line.startswith("tsc stack cost:"):
            tsc_stack_cost = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
            tsc_stack_cost_array.append(tsc_stack_cost / T * 1024)

    # Print the arrays
    print("GC Size Array:", gc_size_array)
    print("Stack Cost Array:", stack_cost_array)
    print("TSC Stack Cost Array:", tsc_stack_cost_array)
    return (
        np.array(gc_size_array),
        np.array(stack_cost_array),
        np.array(tsc_stack_cost_array),
    )


def parse_compute(filename):
    # Initialize arrays
    garbling_time_array = []
    eval_time_array = []

    # Read the file
    with open(filename, "r") as file:
        lines = file.readlines()

    # Process each line
    for line in lines:
        if line.startswith("garbling time:"):
            garbling_time = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
            garbling_time_array.append(garbling_time)
        elif line.startswith("eval time:"):
            eval_time = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
            eval_time_array.append(eval_time)

    # Convert lists to numpy arrays
    garbling_time_array = np.array(garbling_time_array)
    eval_time_array = np.array(eval_time_array)
    print("Garbling Time:", garbling_time_array)
    print("Eval Time:", eval_time_array)

    return garbling_time_array, eval_time_array
