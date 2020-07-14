import numpy as np


def MakeModel():
    S[c, t+1] = np.exp(-l[c, t]) * S[c, t]
    E[c, 1, t+1] = (1 - np.exp(-l[c, t])) * S[c, t]
    for i in range(2, 10):
        E[c, i, t+1] = (1 - p[i-1]) * (1-q[i-1]) * E[c, i-1, t]
    I[c, 1, t+1] = np.sum(np.multiply(np.multiply(p, q), E[c, :, t]))
    I[c, 2, t+1] = (1 - h[c, 1]) * I[c, 1, t]
    I[c, 3, t+1] = (1-h[c, 2]) * I[c, 2, t] + \
        (1 - r) * (1 - h[c, 3]) * I[c, 3, t]
    for j in range(4, 5):
        I[c, j, t+1] = r * (1 - h[c, j-1]) * I[c, j-1, t] + \
            (1 - r) * (1 - h[c, j]) * I[c, j, t]
    R[c, t+1] = R[c, t] + r * I[c, 5, t] + r * I[m, 5, t]
    S[h, t+1] = np.exp(-l[h, t]) * S[h, t]
    E[h, 1, t+1] = (1 - np.exp(-l[h, t])) * S[h, t]
    for i in range(2, 10):
        E[h, i, t+1] = (1 - p[i-1]) * (1-q[i-1]) * E[h, i-1, t]
    I[h, 1, t+1] = np.sum(np.multiply(np.multiply(p, q), E[h, :, t]))
    I[h, 2, t+1] = (1 - h[h, 1]) * I[h, 1, t]
    I[h, 3, t+1] = (1-h[h, 2]) * I[h, 2, t] + \
        (1 - r) * (1 - h[h, 3]) * I[h, 3, t]
    for j in range(4, 5):
        I[h, j, t+1] = r * (1 - h[h, j-1]) * I[h, j-1, t] + \
            (1 - r) * (1 - h[h, j]) * I[h, j, t]
    R[h, t+1] = R[h, t] + r * I[h, 5, t] + r * I[m, 5, t]
    for i in range(2, 10):
        E[m, i, t+1] = (1 - p[i - 1]) * np.floor(q[i - 1] * E)
    I[m, 1, t+1]
    I[m, 2, t+1]
    I[m, 3, t+1]
    for j in range(4, 5):
        I[m, j, t+1]
