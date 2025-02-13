# plot the cost of each component
from basic_cost import *
from sim_cost import (
    cost_recursive,
    time_g_and,
    time_e_and,
    cpu_freq,
    comm_and,
    comm_join,
    time_g_buffer,
    time_e_buffer,
)
import matplotlib.pyplot as plt
from parse_raw import parse_amortized_cost, parse_compute

save_plot = True

w = 64
T = [2**i for i in range(10, 25)]

plt.figure(figsize=(4, 4))

# from nanogram's simulator
nanogram = np.array(
    [
        26784883.0,
        41147730.0,
        52213399.0,
        70949099.0,
        86272094.0,
        112710067.0,
        137789392.0,
        174707282.0,
        202730838.0,
        250616114.0,
        286213645.0,
        350803724.0,
        396293071.0,
        470804247.0,
        526104771.0,
        616127709.0,
    ]
) / (8 * 1024 * 1024)

nanogram_boolean = np.array(
    [
        21532309.0,
        31582984.0,
        39292466.0,
        52140508.0,
        62970761.0,
        80793786.0,
        99679755.0,
        124698495.0,
        144507991.0,
        176653024.0,
        201695783.0,
        246147646.0,
        278421531.0,
        327839644.0,
        366944043.0,
        426362198.0,
    ]
) / (8 * 1024 * 1024)

picogram, picogram_compact, _ = parse_amortized_cost("simd_raw.txt")

pico_no_simd, pico_no_simd_compact, _ = parse_amortized_cost("tsc_raw.txt")

tsc_with_stack_cost = [cost_recursive(t, t, w, 2, "TSC") for t in T]

tsc_with_stack = np.array([cost.comm_MB() for cost in tsc_with_stack_cost])

tsc_stack = np.array([cost.get_stack_cost().comm_MB() for cost in tsc_with_stack_cost])

print(f"For T ranging from 2^10 to 2^24")
print("improvement over nanogram: ", nanogram[: len(T)] / picogram)
print("improvement over tsc", tsc_with_stack / picogram)

width = 64
linear = np.array([2 * 197 * width * t / (8 * 1024 * 1024) for t in T])

plt.plot(np.log2(T), linear, marker="*", label="Linear scan")

plt.plot(
    np.log2(T)[: len(tsc_with_stack)],
    tsc_with_stack,
    marker="o",
    label="TSC",
)

plt.plot(
    np.log2(T),
    nanogram[: len(T)],
    marker="s",
    label="NanoGRAM",
    color="plum",
)

plt.plot(
    np.log2(T)[: len(pico_no_simd)],
    pico_no_simd,
    marker="x",
    label="Pico. no-SIMD",
)


plt.plot(
    np.log2(T),
    picogram,
    marker="d",
    label="PicoGRAM",
    color="red",
)

plt.xlabel("log($N$)", fontsize=12)
plt.ylabel("Cost per Access (MB)", fontsize=12)
plt.legend()

plt.xlim(10, 24)
plt.ylim(0, 100)

if save_plot:
    plt.savefig("cost.pdf", bbox_inches="tight")

# plot a bar chart of the cost breakdown for RAM of size 2^16
plt.figure(figsize=(4, 4))
picogram_boolean = picogram - picogram_compact
pico_no_simd_boolean = pico_no_simd - pico_no_simd_compact
tsc_boolean = tsc_with_stack - tsc_stack
nanogram_stack = nanogram - nanogram_boolean

index = 16 - 10
barWidth = 0.3

picogram_sim_cost = cost_recursive(T[index], T[index], w, 4, "picogram")
# each translation translates to additional 128 bits of comm in the interactive model
interact_comm = picogram_boolean[index] + picogram_sim_cost.translation_g * comm_join / 8 / 1024 / 1024

print(f"At N = 2^{index}")

print("tsc total / interactive comm: ", tsc_with_stack[index] / interact_comm)

print("nanogram total / interactive comm: ", nanogram[index] / interact_comm)

print("picogram total / interactive comm: ", picogram[index] / interact_comm)

plt.bar(
    np.arange(4) - barWidth / 2,
    [
        tsc_stack[index],
        nanogram_stack[index],
        pico_no_simd_compact[index],
        picogram_compact[index],
    ],
    width=barWidth,
    label="Stack",
)

compact_improvement_for_stack = tsc_stack[index] / pico_no_simd_compact[index]
simd_improvement_for_stack = pico_no_simd_compact[index] / picogram_compact[index]
print("compaction improvement for stack: ", compact_improvement_for_stack)
print("simd extra improvement for stack: ", simd_improvement_for_stack)
print(
    "total improvement for stack: ",
    compact_improvement_for_stack * simd_improvement_for_stack,
)

plt.bar(
    np.arange(4) + barWidth / 2,
    [
        tsc_boolean[index],
        nanogram_boolean[index],
        pico_no_simd_boolean[index],
        picogram_boolean[index],
    ],
    width=barWidth,
    label="Other",
)

# name the bars
plt.xticks(np.arange(4), ["TSC", "NanoGRAM", "PicoGRAM\nno-SIMD", "PicoGRAM"])


plt.ylabel("Cost per Access (MB)", fontsize=12)
plt.ylim(0)
plt.legend()
# ax[2].ylim(0, 105)

if save_plot:
    plt.savefig("breakdown.pdf", bbox_inches="tight")


T = [2**i for i in range(10, 17)]

plt.figure(figsize=(4, 4))

# plot the garbler and evaluator's runtime
garble_time, eval_time = parse_compute("simd_compute.txt")

plt.plot(np.log2(T), garble_time / T, marker="o", label="Garbler")

plt.plot(np.log2(T), eval_time / T, marker="s", label="Evaluator")

plt.xlabel("log($N$)", fontsize=12)
plt.ylabel("Time per Access (ms)", fontsize=12)

plt.legend()
plt.ylim(0)

# save the plot
if save_plot:
    plt.savefig("compute.pdf", bbox_inches="tight")

# create a plot that simulate the total runtime for N = 2^15 for each scheme at different bandwidths

# bandwidths = np.array([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
# generate the log scale bandwidths with many points
bandwidths = np.logspace(1, 3.6, 100)


def total_runtime(
    T,
    avg_cost,
    bandwidth,
    total_garble_time_ms=0,
    total_eval_time_ms=0,
    interaction_time_ms=0,
):
    total_comm_time = avg_cost * 8 * 2**20 * 1000 / (1e6 * bandwidth)
    if interaction_time_ms:
        total_eval_time_ms = max(total_eval_time_ms, interaction_time_ms)
    # garbler pipelined with communication
    return max(total_comm_time, total_garble_time_ms / T) + total_eval_time_ms / T


logN = 16
index = logN - 10
garble_thread = 8


def boolean_comm_to_total_eval_ms(comm_MB, T):
    return comm_MB * 8 * 2**20 / comm_and * time_e_and * T / 1e3


def boolean_comm_to_total_garble_ms(comm_MB, T):
    return comm_MB * 8 * 2**20 / comm_and * time_g_and * T / 1e3 / garble_thread


def stack_comm_to_total_eval_ms(comm_MB, T):
    # twice number of buffer than join
    # 3/8 are evaluated
    return comm_MB * 8 * 2**20 / comm_join * 2 * 3 / 8 * time_e_buffer * T / 1e3


def stack_comm_to_total_garble_ms(comm_MB, T):
    # twice number of buffer than join
    return comm_MB * 8 * 2**20 / comm_join * 2 * time_g_buffer * T / 1e3 / garble_thread


linear_garble_time = boolean_comm_to_total_garble_ms(linear[index], T[index])
linear_eval_time = boolean_comm_to_total_eval_ms(linear[index], T[index])

linear_runtime = [
    total_runtime(T[index], linear[index], b, linear_garble_time, linear_eval_time)
    for b in bandwidths
]


def cycles_to_total_ms(cycles, freq, T):
    return cycles / freq * 1e6 * T


tsc_garble_time = (
    cycles_to_total_ms(tsc_with_stack_cost[index].gen_cycles(), cpu_freq, T[index])
    / garble_thread
)  # assume perfect multi-threading

tsc_eval_time = cycles_to_total_ms(
    tsc_with_stack_cost[index].eval_cycles(), cpu_freq, T[index]
)

tsc_runtime = [
    total_runtime(T[index], tsc_with_stack[index], b, tsc_garble_time, tsc_eval_time)
    for b in bandwidths
]

nanogram_garble_time = boolean_comm_to_total_garble_ms(
    nanogram_boolean[index], T[index]
) + stack_comm_to_total_garble_ms(nanogram_stack[index], T[index])

# todo, add cost from stack
nanogram_eval_time = boolean_comm_to_total_eval_ms(
    nanogram_boolean[index], T[index]
) + stack_comm_to_total_eval_ms(nanogram_stack[index], T[index])

nanogram_runtime = [
    total_runtime(
        T[index], nanogram[index], b, nanogram_garble_time, nanogram_eval_time
    )
    for b in bandwidths
]
pico_no_simd_compact_runtime = [
    total_runtime(T[index], pico_no_simd[index], b) for b in bandwidths
]
picogram_runtime = [
    total_runtime(
        T[index],
        picogram[index],
        b,
        garble_time[index],
        eval_time[index],
    )
    for b in bandwidths
]

rtt = 100
interact_round = T[index] * (np.ceil(np.log2(T[index] / 1024)) - 1)
circuit_oram_interactive_runtime = [
    total_runtime(
        T[index],
        interact_comm,
        b,
        0,
        0,
        interact_round * rtt,
    )
    for b in bandwidths
]

bw = 300
picogram_runtime_at_bw = total_runtime(
    T[index], picogram[index], bw, garble_time[index], eval_time[index]
)
nanogram_runtime_at_bw = total_runtime(T[index], nanogram[index], bw, nanogram_garble_time, nanogram_eval_time)
tsc_runtime_at_bw = total_runtime(T[index], tsc_with_stack[index], bw, tsc_garble_time, tsc_eval_time)
interactive_runtime_at_bw = total_runtime(T[index], picogram_boolean[index], bw, 0, 0, interact_round * rtt)

circuit_oram_interactive_runtime_at_bw = total_runtime(
    T[index],
    interact_comm,
    bw,
    0,
    0,
    T[index] * (np.ceil(np.log2(T[index] / 1024)) - 1) * rtt,
)


print("speedup over nanogram: ", nanogram_runtime_at_bw / picogram_runtime_at_bw)

print("speedup over tsc: ", tsc_runtime_at_bw / picogram_runtime_at_bw)

print("speedup over interactive: ", interactive_runtime_at_bw / picogram_runtime_at_bw)

adaptive_garble_time, adaptive_eval_time = parse_compute("simd65536_adaptive.txt")
adaptive_comm, _, _ = parse_amortized_cost("simd65536_adaptive.txt")

picogram_adaptive_runtime = [
    min(
        [
            total_runtime(
                T[index],
                adaptive_comm[i],
                b,
                adaptive_garble_time[i],
                adaptive_eval_time[i],
            )
            for i in range(len(adaptive_comm))
        ]
    )
    for b in bandwidths
]

plt.figure(figsize=(4, 4))

plt.plot(bandwidths, linear_runtime, label="Linear scan")
plt.plot(bandwidths, circuit_oram_interactive_runtime, label="Interactive")
plt.plot(bandwidths, tsc_runtime, label="TSC")
plt.plot(bandwidths, nanogram_runtime, label="NanoGRAM", color="plum")
# plt.plot(bandwidths, pico_no_simd_compact_runtime, label="TSC + Compact")
plt.plot(bandwidths, picogram_runtime, label="PicoGRAM", color="red")
# dashed line for picogram adaptive
plt.plot(
    bandwidths,
    picogram_adaptive_runtime,
    label="Pico. adaptive",
    linestyle=":",
    color="purple",
)

plt.xscale("log")
plt.yscale("log")
plt.xlim(20, 3000)
plt.ylim(50, 30000)

plt.xlabel("Bandwidth (Mbps)", fontsize=12)
plt.ylabel("Runtime per access (ms)", fontsize=12)

# legend on the top right
plt.legend(loc="upper right")

# breakeven bw of picogram and nanogram, tsc, interactive 
for i in range(len(bandwidths)):
    if circuit_oram_interactive_runtime[i] > picogram_runtime[i]:
        print("Breakeven BW of interactive and picogram: ", bandwidths[i])
        break

for i in range(len(bandwidths)):
    if nanogram_runtime[i] < picogram_runtime[i]:
        print("Breakeven BW of nanogram and picogram: ", bandwidths[i])
        break

for i in range(len(bandwidths)):
    if tsc_runtime[i] < picogram_runtime[i]:
        print("Breakeven BW of tsc and picogram: ", bandwidths[i])
        break

if save_plot:
    plt.savefig("runtime.pdf", bbox_inches="tight")
