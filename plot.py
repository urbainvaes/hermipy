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

def plot_projections(series, quad, factor, degree):
    degrees = np.arange(degree + 1)
    series_x = series.project('x')
    series_y = series.project('y')
    quad_x = quad.project('x')
    quad_y = quad.project('y')
    fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)
    quad_x.plot(series_x, degree, factor, ax11)
    ax11.set_title("Marginal distribution in x direction")
    ax12.bar(degrees, series_x.coeffs)
    ax12.set_title("Hermite coefficients of x marginal")
    quad_y.plot(series_y, degree, factor, ax21)
    ax21.set_title("Marginal distribution in y direction")
    ax22.bar(degrees, series_y.coeffs)
    ax22.set_title("Hermite coefficients of y marginal")
    plt.show()
