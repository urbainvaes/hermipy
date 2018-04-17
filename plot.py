import numpy as np

from libhermite import hermite as hm

global quad

def plot_with_quad(series, factor, quad, degree):

    x_visu = quad_visu_x.discretize('x')
    y_visu = quad_visu_y.discretize('x')

    nnodes = [len(nodes_dir) for nodes_dir in quad.nodes]

    solution = quad.eval(series, degree)*factor
    solution_visu_xy = solution.reshape(*nnodes).T

    cont = plt.contourf(x_visu, y_visu, solution_visu_xy, 100)
    plt.title("Eigenvalue: " + str(eigen_value) + ", Minimum: " + str(np.min(solution_visu_xy)))
    plt.colorbar(cont)
    plt.show()
