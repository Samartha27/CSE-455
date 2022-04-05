#include <math.h>
#include <stdlib.h>
#include "image.h"
#include "matrix.h"

// Run an activation function on each element in a matrix,
// modifies the matrix in place
// matrix m: Input to activation function
// ACTIVATION a: function to run
void activate_matrix(matrix m, ACTIVATION a)
{
    int i, j;
    for(i = 0; i < m.rows; ++i){
        double sum = 0;
        for(j = 0; j < m.cols; ++j){
            double x = m.data[i][j];
            if(a == LOGISTIC){
                // TODO
                m.data[i][j] = 1/(1 + expf(-x));
            } else if (a == RELU){
                // TODO
                m.data[i][j] = MAX(0.0,x);
            } else if (a == LRELU){
                // TODO
                m.data[i][j] = MAX(0.1*x,x);
            } else if (a == SOFTMAX){
                // TODO
                m.data[i][j] = expf(x);
            }
            sum += m.data[i][j];
        }
        if (a == SOFTMAX) {
            // TODO: have to normalize by sum if we are using SOFTMAX
            for (j = 0;j<m.cols;j++){
                m.data[i][j] = m.data[i][j]/sum;
            }
        }
    }
}

// Calculates the gradient of an activation function and multiplies it into
// the delta for a layer
// matrix m: an activated layer output
// ACTIVATION a: activation function for a layer
// matrix d: delta before activation gradient
void gradient_matrix(matrix m, ACTIVATION a, matrix d)
{
    int i, j;
    float grad;
    for(i = 0; i < m.rows; ++i){
        for(j = 0; j < m.cols; ++j){
            double x = m.data[i][j];
            // TODO: multiply the correct element of d by the gradient
            if(a == LOGISTIC){
                // TODO
                grad = x*(1-x);
                d.data[i][j] = d.data[i][j] * grad;
                
            } else if (a == RELU){
                // TODO
                if (x>0) {
                    grad = 1;
                    d.data[i][j] = d.data[i][j]*grad;
                }
                else{
                    grad = 0;
                    d.data[i][j] = d.data[i][j]*grad;
                }
                
            } else if (a == LRELU){
                // TODO
                if (x>0) {
                    grad = 1;
                    d.data[i][j] = d.data[i][j]*grad;
                }
                else{
                    grad = 0.1;
                    d.data[i][j] = d.data[i][j]*grad;
                }
            } else if (a == SOFTMAX){
                // TODO
                grad = 1;
                d.data[i][j] = d.data[i][j]*grad;
            }
        }
    }
}

// Forward propagate information through a layer
// layer *l: pointer to the layer
// matrix in: input to layer
// returns: matrix that is output of the layer
matrix forward_layer(layer *l, matrix in)
{

    l->in = in;  // Save the input for backpropagation


    // TODO: fix this! multiply input by weights and apply activation function.
    matrix out = make_matrix(in.rows, l->w.cols);
    out = matrix_mult_matrix(l->in,l->w);
    activate_matrix(out,l->activation); 


    free_matrix(l->out);// free the old output
    l->out = out;       // Save the current output for gradient calculation
    return out;
}

// Backward propagate derivatives through a layer
// layer *l: pointer to the layer
// matrix delta: partial derivative of loss w.r.t. output of layer
// returns: matrix, partial derivative of loss w.r.t. input to layer
matrix backward_layer(layer *l, matrix delta)
{
    // 1.4.1
    // delta is dL/dy
    // TODO: modify it in place to be dL/d(xw)
    gradient_matrix(l->out,l->activation,delta);


    // 1.4.2
    // TODO: then calculate dL/dw and save it in l->dw
    free_matrix(l->dw);
    matrix dw = make_matrix(l->w.rows, l->w.cols); // replace this

    matrix xt = transpose_matrix(l->in);
    dw = matrix_mult_matrix(xt,delta);
    l->dw = dw;
    free_matrix(xt);

    
    // 1.4.3
    // TODO: finally, calculate dL/dx and return it.
    matrix dx = make_matrix(l->in.rows, l->in.cols); // replace this
    matrix wt = transpose_matrix(l->w);
    dx = matrix_mult_matrix(delta,wt);
    free_matrix(wt);


    return dx;
}

// Update the weights at layer l
// layer *l: pointer to the layer
// double rate: learning rate
// double momentum: amount of momentum to use
// double decay: value for weight decay
void update_layer(layer *l, double rate, double momentum, double decay)
{
    // TODO:
    // Calculate Δw_t = dL/dw_t - λw_t + mΔw_{t-1}
    // save it to l->v
    // rate = 0.0001 - 0.01 , lambda = 0.0005 , m = 0.9
    
    matrix wt_temp = axpy_matrix(-decay,l->w,l->dw);
    matrix delta_wt = axpy_matrix(momentum,l->v,wt_temp);
    free_matrix(l->v);
    l->v = delta_wt;
    free_matrix(wt_temp);


    // Update l->w
    l->w = axpy_matrix(rate,l->v,l->w);



    // Remember to free any intermediate results to avoid memory leaks

}

// Make a new layer for our model
// int input: number of inputs to the layer
// int output: number of outputs from the layer
// ACTIVATION activation: the activation function to use
layer make_layer(int input, int output, ACTIVATION activation)
{
    layer l;
    l.in  = make_matrix(1,1);
    l.out = make_matrix(1,1);
    l.w   = random_matrix(input, output, sqrt(2./input));
    l.v   = make_matrix(input, output);
    l.dw  = make_matrix(input, output);
    l.activation = activation;
    return l;
}

// Run a model on input X
// model m: model to run
// matrix X: input to model
// returns: result matrix
matrix forward_model(model m, matrix X)
{
    int i;
    for(i = 0; i < m.n; ++i){
        X = forward_layer(m.layers + i, X);
    }
    return X;
}

// Run a model backward given gradient dL
// model m: model to run
// matrix dL: partial derivative of loss w.r.t. model output dL/dy
void backward_model(model m, matrix dL)
{
    matrix d = copy_matrix(dL);
    int i;
    for(i = m.n-1; i >= 0; --i){
        matrix prev = backward_layer(m.layers + i, d);
        free_matrix(d);
        d = prev;
    }
    free_matrix(d);
}

// Update the model weights
// model m: model to update
// double rate: learning rate
// double momentum: amount of momentum to use
// double decay: value for weight decay
void update_model(model m, double rate, double momentum, double decay)
{
    int i;
    for(i = 0; i < m.n; ++i){
        update_layer(m.layers + i, rate, momentum, decay);
    }
} 

// Find the index of the maximum element in an array
// double *a: array
// int n: size of a, |a|
// returns: index of maximum element
int max_index(double *a, int n)
{
    if(n <= 0) return -1;
    int i;
    int max_i = 0;
    double max = a[0];
    for (i = 1; i < n; ++i) {
        if (a[i] > max){
            max = a[i];
            max_i = i;
        }
    }
    return max_i;
}

// Calculate the accuracy of a model on some data d
// model m: model to run
// data d: data to run on
// returns: accuracy, number correct / total
double accuracy_model(model m, data d)
{
    matrix p = forward_model(m, d.X);
    int i;
    int correct = 0;
    for(i = 0; i < d.y.rows; ++i){
        if(max_index(d.y.data[i], d.y.cols) == max_index(p.data[i], p.cols)) ++correct;
    }
    return (double)correct / d.y.rows;
}

// Calculate the cross-entropy loss for a set of predictions
// matrix y: the correct values
// matrix p: the predictions
// returns: average cross-entropy loss over data points, 1/n Σ(-ylog(p))
double cross_entropy_loss(matrix y, matrix p)
{
    int i, j;
    double sum = 0;
    for(i = 0; i < y.rows; ++i){
        for(j = 0; j < y.cols; ++j){
            sum += -y.data[i][j]*log(p.data[i][j]);
        }
    }
    return sum/y.rows;
}


// Train a model on a dataset using SGD
// model m: model to train
// data d: dataset to train on
// int batch: batch size for SGD
// int iters: number of iterations of SGD to run (i.e. how many batches)
// double rate: learning rate
// double momentum: momentum
// double decay: weight decay
void train_model(model m, data d, int batch, int iters, double rate, double momentum, double decay)
{
    int e;
    for(e = 0; e < iters; ++e){
        data b = random_batch(d, batch);
        matrix p = forward_model(m, b.X);
        fprintf(stderr, "%06d: Loss: %f\n", e, cross_entropy_loss(b.y, p));
        matrix dL = axpy_matrix(-1, p, b.y); // partial derivative of loss dL/dy
        backward_model(m, dL);
        update_model(m, rate/batch, momentum, decay);
        free_matrix(dL);
        free_data(b);
    }
}


// Questions 
//
// 5.2.2.1 Why might we be interested in both training accuracy and testing accuracy? What do these two numbers tell us about our current model?
// TODO
//The training accuracy tells us how well the model performs on our training dataset and whether the model is able 
// to learn complexity of different features. The goal is to reduce
// the training error so that the model learns to predict with higher accuracy. The testing accuracy 
// tells us how well the model performs on the testing dataset which is a generalization of how
// the model performs on data it has never seen before.

// 5.2.2.2 Try varying the model parameter for learning rate to different powers of 10 (i.e. 10^1, 10^0, 10^-1, 10^-2, 10^-3) and
// training the model. What patterns do you see and how does the choice of learning rate affect both the loss during training and 
// the final model accuracy?
// TODO
//              Learning Rate	Iterations	Training Accuracy	Test Accuracy	Loss
//              10	            1000	    Nan                 Nan	            Nan
//              1	            1000	    0.8774              0.8718	        0.449439
//              0.1	            1000	    0.918               0.916           0.2649
//              0.01	        1000	    0.902	            0.907           0.3304
//              0.001	        1000	    0.8585	            0.8689          0.5966
// As the learning rate increases, the model starts learning faster and quickly converges. But the decrease in accuracy
// when the learning rate approaches 1 becomes prominent and is attributed to slower convergence and non convergence when
// the learning rate goes to 10. 



// 5.2.2.3 Try varying the parameter for weight decay to different powers of 10: (10^0, 10^-1, 10^-2, 10^-3, 10^-4, 10^-5).
// How does weight decay affect the final model training and test accuracy?
// TODO
//
//              Decay	Training Accuracy   Test Accuracy       Final Loss
//              1	    0.87421	            0.8803	            0.508432
//              0.1	    0.88635	            0.8915	            0.339479
//              0.01	0.8874	            0.8933	            0.324052
//              0.001	0.8876	            0.8934	            0.322535
//              0.0001	0.8876	            0.8934	            0.322384
//              0	    0.8876	            0.8934	            0.322367

// It can be seen that decreasing the value of decay rate increases the training and the test accuracy.
// The change in the accuracies is very minimal. When the decay increases the training and the test accuracies decreases.
// The model will not be able to effectively learn complexity and hence when decay increases, the accuracies are decreasing.

// 5.2.3.1 Currently the model uses a logistic activation for the first layer. 
//Try using a the different activation functions we programmed. How well do they perform? What's best?
// TODO
// Activation	Loss	    Training_accuracy	Test_accuracy
//LOGISTICS	    0.330491	0.90235	            0.9077
//RELU	        0.132339	0.92315	            0.924
//LRELU	        0.133097	0.920716667	        0.9217
//SOFTMAX	    -nan	    0.098716667	        0.098
// It can be seen that the accuracy of the model varies in the direction  
// SOFTMAX  ->  LOGISTICS  ->  LRELU  ->  RELU(from worst to best). The difference between LRELU and RELU is not very much.
// But the best model for the dataset is RELU. 
//The worst activation function to be used for the input and hidden layers is Softmax since its accuracy is almost 0.
//


// 5.2.3.2 Using the same activation, find the best (power of 10) learning rate for your model. 
//What is the training accuracy and testing accuracy?
// TODO
// Learning rate	Loss	    Training_accuracy	Test_accuracy
// 1	            -nan	    0.098716667	        0.098
// 0.1	            0.070356	0.96115	            0.9564
// 0.01	            0.132339	0.92315	            0.924
// 0.001	        0.41079	    0.864433333	        0.8682
// 0.0001	        1.879566	0.592433333	        0.5915
// 
// It can be seen that the best learning rate is 0.1 (10^(-1)).
// The model has a training accuracy of 0.96115 and testing accuracy of 0.9564
// 


// 5.2.3.3 Right now the regularization parameter `decay` is set to 0.
// Try adding some decay to your model. What happens, does it help? Why or why not may this be?
// TODO
// Decay	Loss	Training_accuracy	Test_accuracy
// 1	    0.178535	0.924516667	    0.9272
// 0.1	    0.076561	0.95685	        0.952
// 0.01	    0.083339	0.961983333	    0.9563
// 0.001	0.092674	0.961933333	    0.9566
// 0.0001	0.064849	0.96245	        0.9586

//It can be seen that adding the decay parameter helps the model. 
//If we choose the decay rate appropriately it can help the model to improve its accuracy.
// In this particular variation it can be seen that reducing the value of the decay improves the accuracy as
// the value of the regularization parameter approaches zero. 
//Having a higher weight decay parameter penalizes the weights more and this reduces overfitting. 
//If we are able to choose the approriate decay along with the appropriate learning rate, 
//the model wil learn to generalize well and obtain better accuracies.
//
//
// 5.2.3.4 Modify your model so it has 3 layers instead of two. 
//The layers should be `inputs -> 64`, `64 -> 32`, and `32 -> outputs`. 
//Also modify your model to train for 3000 iterations instead of 1000. 
//Look at the training and testing error for different values of decay (powers of 10, 10^-4 -> 10^0). Which is best? Why?
// TODO
//
//              Decay	    Training accuracy	Test accuracy	    Loss
//              1	        0.9431	            0.9439	            0.193929
//              0.1	        0.9684	            0.9636	            0.086949
//              0.01	    0.9708	            0.9647	            0.075409
//              0.001	    0.9711	            0.9646	            0.077661
//              0.0001	    0.9709	            0.9644	            0.074904
//              0	        0.9711	            0.9644	            0.075054
//
//  Here,as the decay increases,the training and test accuracies decreases.
// This is because the complexity of the model is decreasng but it should have been increasing to learn the model well.
// Since the model is not overfitting,decreasing the model complexity will decrease the test accuracy.
// Hence the training and test accuracies both are greater for 0 decay.



// 5.3.2.1 How well does your network perform on the CIFAR dataset?
// TODO 
// With a learning rate of 0.01 and decay of 0,
//The model did not perform well and itgave a train accuracy of 0.4558 and a test accuracy of 0.4482.
// After 2000 iterations the error was 1.45.
//



