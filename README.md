# GCN-Denoiser: Mesh Denoising with Graph Convolutional Networks

The project was created in order to provide an opportunity to generate noise for 3D meshes. This is necessary to create a dataset for learning-based denoising algorithms. This project provides a comprehensive library with 3D mesh noise generation methods.

## Code

### Prerequisites:

- Hardware: Personal computer;
- Environments: Windows system.

### Third Party Library:

- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page);
- [Glm](https://github.com/g-truc/glm);
- [Lz4](https://github.com/lz4/lz4);
- [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/);
- [PyBind11](https://github.com/pybind/pybind11).

### Denoising Interface:

Executable demo, the corresponding code, and some sampled meshes are supplied.

- For .pyd testing the Jupyter notebook is required with python 3.8. If you want to use another python version you need to rebuild the library and change the pyhton version in compiler include libraries and linker additional dependencies:

- For code, Visual Studio 2019 is required.

## Algorithms provided

The algorithms implemented are divided into two groups: node noise algorithms and topology noise algorithms.

- Each function reads an obj file from specified directory and writes another obj file to the directory where original obj file is placed;
- Each function has a seed argument. Default value is 0;
- Each function first argument is modelPath which specifies an absolute obj path. It is required to specify this argument for proper function working.

Each node noise generation function has a noiseDirection argument of int type which specifies the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Node normal is calculated as a mean of adjacent faces normals.

* `GenerateGaussianNoisePy(string modelPath, double sigma, int noiseDirection, int seed)` -- generates node noise distributed by Gaussian PDF: $P(x | \sigma) = \frac{1}{\sigma \sqrt{2 \pi}} e^{-\frac{1}{2} \frac{x^{2}}{\sigma^{2}}}$.

  - string modelPath -- absolute .obj path. Must be specified;
  - double sigma -- standard deviation of the Gaussian noise model. Defaultvalue is 0.2;
  - int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateImpulsiveNoisePy(string modelPath, double sigma, double verticesPortion, int noiseDirection, int seed)` -- generates node noise distributed by Gaussian PDF. The portion of nodes to be noised is specified.

  - string modelPath -- absolute .obj path. Must be specified;
  - double sigma -- standard deviation of the Gaussian noise model. Defaultvalue is 0.2;
  - verticesPortion -- a portion of vertices to be noised. Default value is 0.2. Maximum value is 1;
  - int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateExponentialNoisePy(string meshDir, double lambda, int noiseDirection, int seed)` -- generates node noise distributed by Exponential PDF: $P(x | \lambda) = \lambda e^{-\lambda x}$.

  - string modelPath -- absolute .obj path. Must be specified;
  - double lambda -- Exponential distribution $\lambda$ rate parameter. Default value is 7.0;
  - int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateExtremeValueNoisePy(string meshDir, double a, double b, int noiseDirection, int seed)` -- generates node noise distributed by Generalized extreme value PDF: $P(x | a, b) = \frac{1}{b} \exp{\frac{a-x}{b} - \exp{\frac{a-x}{b}}}$.

  - string modelPath -- absolute .obj path. Must be specified;
  - double a -- location parameter. Default value is 0.0;
  - double b -- scale parameter. Default value is 0.3;
  - int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateGammaNoisePy(string meshDir, double alpha, double beta, int noiseDirection, int seed)` -- generates node noise distributed by Gamma PDF: $P(x | \alpha, \beta) = \frac{e^{-x/\beta}}{\beta^{\alpha} \Gamma(\alpha)} \cdot x^{\alpha - 1}$, where $\Gamma(\alpha)$ is Gamma function.

  - string modelPath -- absolute .obj path. Must be specified;
  - double alpha -- shape parameter. Default value is 0.1;
  - double beta -- rate parameter. Default value is 0.3;
  - int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateLaplaceNoisePy(string meshDir, double mu, double b, int noiseDirection, int seed)` -- generates node noise distributed by Laplace PDF: $P(x | \mu, b) = \frac{1}{2b} \exp{-\frac{|x - \mu|}{b}}$.

  - string modelPath -- absolute .obj path. Must be specified;
  - double mu -- location parameter. Default value is 0.0;
  - double b -- scale parameter. Default value is 0.3;
  - int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateLogNormalNoisePy(string meshDir, double m, double s, int noiseDirection, int seed)` -- generates node noise distributed by Lognormal PDF: $P(x | m, s) = \frac{1}{sx \sqrt{2\pi}}\exp{-\frac{(\ln{x} - m)^{2}}{ 2s^{2}}}$.

  - string modelPath -- absolute .obj path. Must be specified;
  - double m -- mean of the natural logarithm of x. Default value is -5.0;
  - double s -- standard deviation of the natural logarithm of x. Default value is 40.0;
  - int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateUniformNoisePy(string meshDir, double a, double b, int noiseDirection, int seed)` -- generates node noise distributed by Uniform PDF: $P(x | a, b) = \frac{1}{b - a}$.

  - string modelPath -- absolute .obj path. Must be specified;
  - double a -- lower bound. Default value is 0.2;
  - double b -- upper bound. Default value is 0.4;
  - int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateWeibullNoisePy(string meshDir, double a, double b, int noiseDirection, int seed)` -- generates node noise distributed by Uniform PDF: $P(x | a, b) = \frac{a}{b} (\frac{x}{b})^{a-1} \exp{-(\frac{x}{b})^{a}}$.

  - string modelPath -- absolute .obj path. Must be specified;
  - double a -- shape parameter. Default value is 1.0;
  - double b -- scale parameter. Default value is 0.2;
  - int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
  - int seed -- random component determination. Default value is 0;

Each topology noise generation function has a meshNoisedDir argument of string type which specifies an absolute output obj path.

* `GenerateRandomVerticesTopologyNoisePy(string meshDir, string meshNoisedDir, double percentage, int seed)` -- Removes random verteces with adjacent faces.

  - string modelPath -- absolute input .obj path. Must be specified;
  - string meshNoisedDir -- absolute output .obj path. Must be specified;
  - double percentage -- percentage of nodes to be removed. Default value is 5.0. Maximum value is 100;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateRandomFacesTopologyNoisePy(string meshDir, string meshNoisedDir, double percentage, int seed)` -- Removes randomly selected faces.

  - string modelPath -- absolute input .obj path. Must be specified;
  - string meshNoisedDir -- absolute output .obj path. Must be specified;
  - double percentage -- percentage of faces to be removed. Default value is 5.0. Maximum value is 100;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateOneRandomClusterTopologyNoisePy(string meshDir, string meshNoisedDir, double standadDeviation, int maxDistance, double divider, int seed)` -- Removes one normally distributed cluster with randomly selected center. The probability of node removing is defined by Gaussian PDF:
  $$
  P(x | \sigma, D, h) = \left\{
  \frac{1}{\sigma \sqrt{2 \pi}} e^{-\frac{1}{2} \frac{x^{2}}{h^{2}\sigma^{2}}} & \mbox{ if x \le D}\\
  0 & \mbox{if x > D}
  $$

  - string modelPath -- absolute input .obj path. Must be specified;
  - string meshNoisedDir -- absolute output .obj path. Must be specified;
  - double standadDeviation -- percentage of faces to be removed. Default value is 5.0. Maximum value is 100;
  - int maxDistance -- 
  - double divider --
  - int seed -- random component determination. Default value is 0;
  
* ghjgj
* sdfsdf
* sdfsdf
dfdfd
dfdf
* 


## Examples

```bash
import NoiseGenerator
modelPath = 'absolute_path/model.obj'
```

```bash
NoiseGenerator.GenerateGaussianNoisePy(modelPath, 0.2, 0)
NoiseGenerator.GenerateImpulsiveNoisePy(modelPath, 0.2, 0.5, 0)
NoiseGenerator.GenerateExponentialNoisePy(modelPath, 3, 0)
NoiseGenerator.GenerateExtremeValueNoisePy(modelPath, 0.0, 0.3, 0)
NoiseGenerator.GenerateGammaNoisePy(modelPath, 0.1, 0.2, 0)
NoiseGenerator.GenerateLaplaceNoisePy(modelPath, 0.0, 0.3, 0)
NoiseGenerator.GenerateLogNormalNoisePy(modelPath, -5, 40, 0)
NoiseGenerator.GenerateUniformNoisePy(modelPath, 0.2, 0.4, 0)
NoiseGenerator.GenerateWeibullNoisePy(modelPath, 1, 0.2, 0)
```

```bash
import NoiseGenerator
modelPath = 'absolute_path/model.obj'
modelNoisedPath = 'absolute_path/model_noised.obj'
```

```bash
NoiseGenerator.GenerateRandomVerticesTopologyNoisePy(modelPath, modelNoisedPath, 5)
NoiseGenerator.GenerateOneRandomClusterTopologyNoisePy(modelPath, modelNoisedPath, 0.5, 50, 40.0)
NoiseGenerator.GenerateSetOfRandomClustersTopologyNoisePy(modelPath, modelNoisedPath, 0.5, 50, 40.0, 4, 40, 0)
NoiseGenerator.GenerateSetOfRandomClustersTopologyNoisePy(modelPath, modelNoisedPath, 0.5, 30, 20.0, 5, 40, 1)
NoiseGenerator.GenerateSetPatchTopologyNoisePy(modelPath, modelNoisedPath, 0.005, 0.01, 0.01, 0.01, 0.05, 3, 30, 2, 3, 16)
```



