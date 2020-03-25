import numpy as np
import matplotlib.pyplot as plt
A = np.array([[1,0],[1,1],[1,3],[1,4]])
b = np.array([0,1,2,5])
xi = np.array([1,1])
At = np.transpose(A)


for i in range(25):
    plt.text(xi[0], xi[1], str(i))
    plt.plot(xi[0], xi[1])
    best = b - np.matmul(A,xi)
    grad = np.matmul(At,best)
    xi += grad


plt.show()
