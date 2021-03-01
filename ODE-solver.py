import autograd.numpy as np
from autograd import grad, elementwise_grad
import autograd.numpy.random as npr
from matplotlib import pyplot as plt

#define sigmoid function
def sigmoid(z):
    return 1/(1+np.exp(-z))

#compute and find output of NN
def neural_network(params, x):
    w_hidden = params[0]
    w_output = params[1]

    num_values = np.size(x)
    x=x.reshape(1,num_values)
    #input layer
    x_input = x
    #hidden layer
    x_input = np.concatenate((np.ones((1,num_values)), x_input), axis = 0)

    z_hidden = np.matmul(w_hidden, x_input)
    x_hidden = sigmoid(z_hidden)
    #output layer
    x_hidden = np.concatenate((np.ones((1,num_values)), x_hidden), axis = 0)
    z_output = np.matmul(w_output, x_hidden)

    return z_output

#g_t
def g_trial(x,params, g0 = 10):
    return g0 + x*neural_network(params,x)

def g(x, g_trial, gamma = 2):
    return -gamma*g_trial

def cost_function(P, x):
    g_t = g_trial(x,P)
    #dg/dx
    dg_dx = elementwise_grad(g_trial, 0)(x,P)
    #(dg_dx-g(x))^2
    matrix_squared_error = (dg_dx - g(x,g_t))**2
    mean_squared_error = np.sum(matrix_squared_error)/np.size(matrix_squared_error)

    return mean_squared_error
#function used to optimize P
def optimize_P(x, num_neurons_hidden, num_iter, lmb):
    #Initial random value to Weights and Bias
    p0 = npr.randn(num_neurons_hidden, 2)
    p1 = npr.randn(1,num_neurons_hidden+1)
    P = [p0,p1]

    print('Initial cost',cost_function(P,x))
    cost_function_grad = grad(cost_function,0) 

    for i in range(num_iter):
        cost_grad = cost_function_grad(P,x)
        P[0] = P[0] - lmb*cost_grad[0]
        P[1] = P[1] - lmb*cost_grad[1]

    print('Final cost',cost_function(P,x))
    return P

def g_actual(x, gamma = 2, g0 = 10):
    return g0*np.exp(-gamma*x)
if __name__ == '__main__':
    npr.seed(16)

    training_size = 10
    training_example = np.linspace(0,1,training_size)

    num_hidden_neurons = 10
    num_iter = 20000
    lmb = 0.001

    print('size of input vector:',training_size)
    print('number of hidden neurons:',num_hidden_neurons)
    print('epochs:',num_iter)
    print('alpha:',lmb)

    P = optimize_P(training_example, num_hidden_neurons, num_iter, lmb)

    result_example = np.linspace(0,1.5,100)
    actual_res = g_actual(result_example)
    training_res = g_trial(result_example, P)

    plt.plot(result_example, actual_res)
    plt.plot(result_example, training_res[0,:])
    plt.legend(['actual_result','training_result'])
    plt.show()