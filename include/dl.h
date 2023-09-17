// NeuralNetwork.hpp
#include <iostream>
#include <fstream>      // std::ifstream
#include <vector>
#include <time.h>
#include "Eigen/Eigen"

namespace dl
{
// use typedefs for future ease for changing data types like : float to double
typedef float Scalar;
typedef Eigen::MatrixXf Matrix;
typedef Eigen::RowVectorXf RowVector;
typedef Eigen::VectorXf ColVector;
typedef std::vector<Matrix*> Data;
typedef std::function<Scalar(Scalar)> TAcFun;

// neural network implementation class!
class NeuralNetwork {
public:


	// constructor
	NeuralNetwork();
	NeuralNetwork(std::vector<unsigned> topology, 
		   	    Scalar learningRate = Scalar(0.005),
			    unsigned batch_size = unsigned(100));
	NeuralNetwork(std::vector<Matrix> & init_weights,  
			    std::vector<Matrix> & init_biases, 
			    std::vector<TAcFun> & init_activations);

	// destructor
	~NeuralNetwork();

	//activation function for value z of the neuron_no th neuron in the layer_no th layer.
	Scalar activationFun (Scalar z, unsigned layer_no, unsigned neuron_no = 0);

	Scalar activationDerivativeFun (Scalar z, unsigned layer_no, unsigned neuron_no = 0);

	std::vector<TAcFun> activations;

	// function for forward propagation of data
	void propagateForward(Matrix& input, int f_train);

	// function for backward propagation of errors made by neurons
	void propagateBackward(Matrix& output);

	// function to calculate errors made by neurons in each layer
	void calcErrors(Matrix& output);

	// function to update the weights of connections
	void updateWeights();

	// function to train the neural network give an array of data points
	void train(std::vector<Matrix*> data, std::vector<Matrix*> output_data);

	Scalar cost(Matrix & expected);

	void cost_gradient(Matrix & derivatives, Matrix & expected);

	//cost function
	//void cost ();

	//cost function derivative	
	//dfvoid 	
	// storage objects for working of neural network
	/*
		use pointers when using std::vector<Class> as std::vector<Class> calls destructor of
		Class as soon as it is pushed back! when we use pointers it can't do that, besides
		it also makes our neural network class less heavy!! It would be nice if you can use
		smart pointers instead of usual ones like this
		*/
    	std::vector<unsigned> topology;
    	std::vector<Matrix*> layers;//H_i, which are input, hidden, and output layers.
	std::vector<Matrix*> gradients_H; // stores the gradients of cost C with respect to H_i

	std::vector<Matrix*> dfs; //f'(z) activation Derivatives

	std::vector<Matrix*> weights; // W_i
	std::vector<Matrix*> gradients_W;
	std::vector<Matrix*> gradients_W_batch;

	std::vector<Matrix*> biases; //B_i
	std::vector<Matrix*> gradients_B;
	std::vector<Matrix*> gradients_B_batch;
	Scalar learningRate;
	unsigned batch_size;
};

struct Model1
{
     NeuralNetwork nn; 
     Model1();
     Model1(std::vector<unsigned>);
     Model1(std::vector<Matrix> & init_weights, std::vector<Matrix> & init_biases, std::vector<TAcFun> & init_activations);

     void forward(Matrix & input, int f_train);
	Matrix & getOutput();
	Matrix & getInput();
};
struct Model2
{
     NeuralNetwork nn00; 
     NeuralNetwork nn10;
     NeuralNetwork nn11;
     Model2();
     Model2(std::vector<unsigned>, std::vector<unsigned>, std::vector<unsigned>);
     Matrix & forward(Matrix & input, int f_train);
	Matrix & getOutput(unsigned i = 1);
	Matrix & getInput();
};

//Type 1 autoencoder generative adversarial networks
struct AcGan1
{
	Model1 E; 
	Model1 G1;
	Model1 G2;
	Model2 D1;
	Model2 D2;

	AcGan1(std::vector<unsigned> E_topology,
		  std::vector<unsigned> G_topology,
		  std::vector<unsigned> D_topology1,
		  std::vector<unsigned> D_topology2,
		  std::vector<unsigned> D_topology3
		  );
}; 

//Type 2 autoencoder generative adversarial networks
struct Nn2AncSVParms{
    std::vector<Matrix> E1_weights;
    std::vector<Matrix> G1_weights;
    std::vector<Matrix> G2_weights;
    std::vector<Matrix> D1_weights;
    std::vector<Matrix> D2_weights;

    std::vector<Matrix> E1_biases;
    std::vector<Matrix> G1_biases;
    std::vector<Matrix> G2_biases;
    std::vector<Matrix> D1_biases;
    std::vector<Matrix> D2_biases;

    std::vector<TAcFun> E1_activations;
    std::vector<TAcFun> G1_activations;
    std::vector<TAcFun> G2_activations;
    std::vector<TAcFun> D1_activations;
    std::vector<TAcFun> D2_activations;
    Nn2AncSVParms();
};

struct AcGan2
{
	Model1 E1; 
	Model1 G1;
	Model1 G2;
	Model1 D1;
	Model1 D2;

	AcGan2(std::vector<unsigned> E_topology,
		  std::vector<unsigned> G_topology,
		  std::vector<unsigned> D_topology
		  );
	AcGan2(std::vector<Matrix> & E1_weights, std::vector<Matrix> & E1_biases, std::vector<TAcFun> &,
		  std::vector<Matrix> & G1_weights, std::vector<Matrix> & G1_biases, std::vector<TAcFun> &,
		  std::vector<Matrix> & G2_weights, std::vector<Matrix> & G2_biases, std::vector<TAcFun> &,
		  std::vector<Matrix> & D1_weights, std::vector<Matrix> & D1_biases, std::vector<TAcFun> &,
		  std::vector<Matrix> & D2_weights, std::vector<Matrix> & D2_biases, std::vector<TAcFun> &);
}; 

void ReadCSV(std::string filename, Data & data);
Scalar MSE(NeuralNetwork & nn, std::vector<Matrix*> input_data, std::vector<Matrix*> output_data);

extern AcGan2 nn2_anc_sv;
extern Nn2AncSVParms acgan_parms;
}
