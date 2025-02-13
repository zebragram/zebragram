import numpy as np
import params


def stack_cost(T):
    d = int(np.log2(T))
    n_join = T * (d - 1)
    return n_join * params.k_ddh


def standard_stack_cost(T, w):
    d = int(np.log2(T))
    n_join = T * (d - 1)
    return n_join * params.k_aes * w


def tsc_stack_cost_helper(m, w, T):
    if m <= 0:
        assert T <= 0
        return 0
    return (
        m * 3 + m // 2 * 4 + m * w * 4  # line14  # line 19 and 20  # line 31-33
    ) * (params.k_aes + 1) + tsc_stack_cost_helper(m // 2, w * 2, (T - 3) // 2)


def tsc_stack_cost(T, w):
    return tsc_stack_cost_helper(T, w, T // 2)


def counter_cost(T):
    d = int(np.log2(T))
    n_half_and = (T - 1) * (d - 1)
    return n_half_and * params.k_aes


def stack_with_counter_cost(T):
    return counter_cost(T) + stack_cost(T)


def standard_stack_with_counter_cost(T, w):
    return standard_stack_cost(T, w) + counter_cost(T)


def switch_cost(T, way=2):
    return stack_with_counter_cost(T) * way


def tach_cost(T, w):
    return T * w * (params.k_aes + params.sigma + params.k_ddh * 2)


def switch_with_tach_cost(T, w, way=2):
    return switch_cost(T, way) + tach_cost(T, w)


def standard_switch_cost(T, w, way=2):
    return standard_stack_cost(T, w) * way


def tsc_switch_cost(T, w, way=2):
    return tsc_stack_cost(T, w) * way


def pos_scan_cost(T, addr_width=20):
    d = int(np.log2(T))
    return T * T * (addr_width + d) * (params.k_aes * 1.5 + 5)  # one hot garbling


def bkt_cost(T):
    d = int(np.log2(T))
    tree_cost = 0
    for i in range(d):
        tree_cost += 2**i * switch_cost(T // 2**i)
    return pos_scan_cost(T) + tree_cost


def bkt_cost_with_tach(T, w):
    return bkt_cost(T) + tach_cost(T, w)


def standard_bkt_cost(T, w):
    d = int(np.log2(T))
    tree_cost = 0
    for i in range(d):
        tree_cost += 2**i * standard_switch_cost(T // 2**i, w + d)
    return pos_scan_cost(T) + tree_cost


def linear_scan_cost(T, w):
    addr_width = int(np.log2(T))
    return pos_scan_cost(T, addr_width / 2) * 2 * w  # read / write


if __name__ == "__main__":
    for i in range(1, 30):
        intercept = tsc_stack_cost(2**i, 0)
        perbit = tsc_stack_cost(2**i, 1) - intercept
        print("{" + str(intercept) + ", " + str(perbit) + "},")
