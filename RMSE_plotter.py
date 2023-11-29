import numpy as np
import dataholder as data
import cv2


def solve_antoine(T, component):
    """Takes T (K) and component (1 or 2) and return P_sat (bar)"""
    if component == 1:
        if T <= data.antoine_range_c1_1[1]:
            return pow(10, data.antoine_a_c1_1 - (data.antoine_b_c1_1 / (data.antoine_c_c1_1 + T)))
        return pow(10, data.antoine_a_c1_2 - (data.antoine_b_c1_2 / (data.antoine_c_c1_2 + T)))
    if T <= data.antoine_range_c2_1[1]:
        return pow(10, data.antoine_a_c2_1 - (data.antoine_b_c2_1 / (data.antoine_c_c2_1 + T)))
    return pow(10, data.antoine_a_c2_2 - (data.antoine_b_c2_2 / (data.antoine_c_c2_2 + T)))

def solve_P_25(delta_lambdas, Temp):
    """Takes the Wilson interaction energy parameters and T (C) and returns the sum of the squared errors"""
    R = 8.3145
    T = 273.15 + Temp
    if Temp == 25:
        V1 = 8.16512E-05
        V2 = 1.80542E-05
    else:
        V1 = 8.43256E-05
        V2 = 1.82186E-05

    l12_11, l21_22 = delta_lambdas
    op = 0

    if Temp == 25:
        x_array = data.x1_25c
        p_array = data.P_25c
    else:
        x_array = data.x1_50c
        p_array = data.P_50c
    for i in range(len(x_array)):
        x1 = x_array[i]
        x2 = 1 - x1
        L12 = (V2 / V1) * np.exp(-l12_11 / (R * T))
        L21 = (V1 / V2) * np.exp(-l21_22 / (R * T))
        x1_x2L12 = x1 + x2 * L12
        x2_x1L21 = x2 + x1 * L21
        bracket = (L12 / x1_x2L12) - (L21 / x2_x1L21)
        y1 = np.exp(-np.log(x1_x2L12) + x2 * bracket)
        y2 = np.exp(-np.log(x2_x1L21) - x1 * bracket)

        P_solved = (y1 * x1 * solve_antoine(T, 1) + y2 * x2 * solve_antoine(T, 2)) * 100000
        op += (p_array[i] - P_solved) ** 2
    return op


minval = 10000000
maxval = 0

xrange = 50000
yrange = 50000
xstep = 100
ystep = 100
xmn = -25000
ymn = -25000

temperature = 25

if temperature == 25:
    length = len(data.P_25c)
else:
    length = len(data.P_50c)

img = np.zeros((int(np.ceil(yrange / ystep)), int(np.ceil(xrange / xstep))))
# calculate the RMSE for all parameter pairs and find the minimum and maximum RMSE values
for u in range(0, yrange, ystep):
    for v in range(0, xrange, xstep):
        val = np.sqrt(solve_P_25([v + xmn, u + ymn], temperature) / length)
        if val < minval:
            minval = val
        if val > maxval:
            maxval = val
        img[u // ystep, v // xstep] = val
print(minval)
print(maxval)
print(np.log(maxval - minval))

log_range = np.log(maxval - minval)
img2 = np.zeros((int(np.ceil(yrange / ystep)), int(np.ceil(xrange / xstep)), 3))
# convert the RMSE matrix to a grayscale image
for u in range(0, yrange, ystep):
    for v in range(0, xrange, xstep):
        try:
            col = abs(1 - (np.log((img[u // ystep, v // xstep] - minval)) / log_range)) * 255
            img2[u // ystep, v // xstep] = [col, col, col]
        except:
            print(u, v)
            print(u // ystep, v // xstep)
            print(img[u // ystep, v // xstep] - minval)
            print(np.log(img[u // ystep, v // xstep] - minval))
            raise Exception

# save the image
cv2.imwrite(f"{temperature}c_{xrange//1000}k_output_1.png", img2)
