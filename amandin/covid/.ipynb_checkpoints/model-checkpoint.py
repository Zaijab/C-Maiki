import numpy as np

# Model Described in Paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1691475/


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
    for i in range(1, 10, 1):
        I[m, 1, t+1] = I[m, 1, t+1] + p[i](q[i]*E[c, i, t] + E[m, i, t])
    I[m, 2, t+1] = h[c, 1] * I[c, 1, t] + h[h, 1] * I[h, 1, t] + I[m, 1, t]
    I[m, 3, t+1] = h[c, 2] * I[c, 2, t] + h[h, 2] * I[h, 2, t] + I[m, 2, t] + \
        (1-r) * (h[c, 3] * I[c, 3, t] + h[h, 3] * I[h, 3, t] + I[m, 1, t])
    for j in range(4, 5):
        I[m, j, t+1] = r * (h[c, j-1] * I[c, j-1, t] + h[h, j-1] * I[h, j-1, t] + I[m, j-1, t]
                            ) + (1-r) * (h[c, j] * I[c, j, t] + h[h, j] * I[h, j, t] + I[m, j, t])
    return [S, E, I, R]


def Variables(t):
    S = zeros(1, 2, t+1)
    E = zeros(14, 3, t+1)
    I = zeros(5, 4, t+1)
    R = zeros(1, 2, t+1)

    lambda=zeros(1, 2, t)
    p = zeros(14)
    q = zeros(14)
    h = zeros(5, 2)

    R(1, 1, 1) = 0
    R(1, 2, 1) = 0
    for i = 1:
        14
        E(i, 1, 1) = 3
    end
    E(1, 2, 1) = 0
    E(1, 3, 1) = 0
    I(1, 1, 1) = 1
    I(2, 1, 1) = 0
    I(3, 1, 1) = 0
    I(4, 1, 1) = 0
    I(5, 1, 1) = 0
    I(1, 2, 1) = 0
    I(1, 3, 1) = 0
    I(2, 3, 1) = 0
    I(3, 3, 1) = 0
    I(4, 3, 1) = 0
    I(5, 3, 1) = 0
    I(1, 4, 1) = 0
    S(1, 2, 1) = 15000
    S(1, 1, 1) = 952712-sum(E(: , 1, 1))-sum(I(: , 1, 1))-sum(I(: , 3, 1))-S(1, 2, 1)-...
    R(1, 1, 1)

    p(1) = 0.002*3
    p(2) = 0.005*3
    p(3) = 0.01*3
    p(4) = 0.04*3
    p(5) = 0.08*3
    p(6) = 0.04*3
    p(7) = 0.01*3
    p(8) = 0.005*3
    p(9) = 0.002*3
    p(10) = 0.002*3
    p(11) = 0.002*3
    p(12) = 0.001*3
    p(13) = 0.001*3
    p(14) = 0
    q(1) = 0.0
    q(2) = 0.0
    q(3) = 0.0
    q(4) = 0.0
    q(5) = 0.0
    q(6) = 0.0
    q(7) = 0.0
    q(8) = 0.0
    q(9) = 0.0
    q(10) = 0.0
    q(11) = 0.0
    q(12) = 0.0
    q(13) = 0.0
    q(14) = 0.0
    h(1, 1) = .2
    h(1, 2) = .01
    h(2, 1) = .5
    h(2, 2) = .01
    h(3, 1) = .9
    h(3, 2) = .01
    h(4, 1) = 1
    h(4, 2) = .01
    h(5, 1) = 1
    h(5, 2) = .01
    r = 0.48

    return [S, E, I, R, lambda, p, q, h, r]
