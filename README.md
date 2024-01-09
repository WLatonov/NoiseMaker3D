# Noisemaker 3D: comprehensive framework for mesh noise generation

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

All the third party libraries must be placed into deps folder in CPPNoiseGenerator_PyLib folder.

### Denoising Interface:

Executable demo, the corresponding code, and some sampled meshes are supplied.

- For .pyd testing the Jupyter notebook is required with python 3.8. If you want to use another python version you need to rebuild the library and change the pyhton version in compiler include libraries and linker additional dependencies:

<img src="/PICS/compiler.PNG" style="zoom:60%;" />
<img src="/PICS/linker_1.PNG" style="zoom:60%;" />
<img src="/PICS/linker_1.PNG" style="zoom:60%;" />

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

* `GenerateCauchyNoise(string meshDir, double x0, double gamma, int noiseDirection, int seed)` -- generates node noise distributed by PDF: $(x | x_{0}, \gamma) = \frac{1}{\pi \gamma \[ 1+ \frac{(x - x_{0})^{2}}{\gamma^{2}} \]}$.

  - string modelPath -- absolute .obj path. Must be specified;
  - double x0 -- shape parameter. Default value is 0.0;
  - double gamma -- scale parameter. Default value is 0.5;
  - int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
  - int seed -- random component determination. Default value is 0;

* `GenerateFisherNoise(string meshDir, double d1, double d2, int noiseDirection, int seed)` -- generates node noise distributed by PDF: $(x | d_{1}, d_{2}) = \sqrt{ \frac{(d_{1} x)^{d_{1}} d_{2}^{d_{2}}}{(d_{1}x + d_{2})^{d_{1} + d_{2}} } }/ (x B(\frac{d_{1}}{2}, \frac{d_{2}}{2}))$. Here $B$ is beta function.

  - string modelPath -- absolute .obj path. Must be specified;
  - double d1 -- first degree of freedom parameter. Default value is 1.0;
  - double d2 -- second degree of freedom parameter. Default value is 1.0;
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
  
* `GenerateOneRandomClusterTopologyNoisePy(string meshDir, string meshNoisedDir, double standadDeviation, int maxDistance, double divider, int seed)` -- Removes one normally distributed cluster with randomly selected center. The probability of node removing is defined by Gaussian PDF: $P(x | \sigma, h) = \frac{1}{\sigma \sqrt{2 \pi}} e^{-\frac{1}{2} \frac{x^{2}}{h^{2}\sigma^{2}}}$ if $x \le D$ and $P(x | \sigma, h) = 0$ if $x > D$, where x is a number of edges between node and cluster center.

  - string modelPath -- absolute input .obj path. Must be specified;
  - string meshNoisedDir -- absolute output .obj path. Must be specified;
  - double standadDeviation -- $\sigma$ for probability of node removing. Default value is 0.5;
  - int maxDistance -- $D$, defines maximal radius of cluster. Default value is 50;
  - double divider -- $h$, defunes a scale of Gaussian PDF used for removing probability calculation. Default value is 40.0;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateOneSpecifiedClusterTopologyNoisePy(string meshDir, string meshNoisedDir, double standadDeviation, int maxDistance, double divider, int clusterCenterVertexId, int seed)` -- Removes one normally distributed cluster with randomly selected center. The probability of node removing is defined by Gaussian PDF: $P(x | \sigma, h) = \frac{1}{\sigma \sqrt{2 \pi}} e^{-\frac{1}{2} \frac{x^{2}}{h^{2}\sigma^{2}}}$ if $x \le D$ and $P(x | \sigma, h) = 0$ if $x > D$, where x is a number of edges between node and cluster center.

  - string modelPath -- absolute input .obj path. Must be specified;
  - string meshNoisedDir -- absolute output .obj path. Must be specified;
  - double standadDeviation -- $\sigma$ for probability of node removing. Default value is 0.5;
  - int maxDistance -- $D$, defines maximal radius of cluster. Default value is 50;
  - double divider -- $h$, defunes a scale of Gaussian PDF used for removing probability calculation. Default value is 40.0;
  - int clusterCenterVertexId -- vertex id which is selected to be a cluster center. Default value is 0;
  - int seed -- random component determination. Default value is 0;
  
* `GenerateSetOfRandomClustersTopologyNoisePy(string meshDir, string meshNoisedDir, double standadDeviation, int maxDistance, double divider, int clustersNumber, int interClusterDistance, int makeClustersFar, int seed)` -- Removes set of normally distributed clusters. For each cluster the probability of node removing is defined by Gaussian PDF: $P(x | \sigma, h) = \frac{1}{\sigma \sqrt{2 \pi}} e^{-\frac{1}{2} \frac{x^{2}}{h^{2}\sigma^{2}}}$ if $x \le D$ and $P(x | \sigma, h) = 0$ if $x > D$, where x is a number of edges between node and cluster center.

  - string modelPath -- absolute input .obj path. Must be specified;
  - string meshNoisedDir -- absolute output .obj path. Must be specified;
  - double standadDeviation -- $\sigma$ for probability of node removing. Default value is 0.5;
  - int maxDistance -- $D$, defines maximal radius of cluster. Default value is 50;
  - double divider -- $h$, defunes a scale of Gaussian PDF used for removing probability calculation. Default value is 40.0;
  - int clustersNumber -- the number of clusters to be removed. Default value is 3;
  - int interClusterDistance -- the minimal distance between clusters centers. Default value is 100;
  - int makeClustersFar -- defines if distance between clusters centers must be maximal (0 -- minimize distance between centers, 1 -- maximize distance between centers). Default value is 0;
  - int seed -- random component determination. Default value is 0;

* `GenerateSetPatchTopologyNoisePy(string meshDir, string meshNoisedDir, double c1, double c2, double c3, double c4, double c5, int type, double percentage, int ringsNumber, int ringsNumberToDelete, int regionRadius, int seed)` -- Removes specified patches selected by restricting constants. The special patches are calculated for each faces in mesh. The patches are classified into four types: Plane, Edge, Corner and Transitional. Any types of patches can be selected to delete. The eigen values $(\lambda_{1}, \lambda_{2}, \lambda_{3})$ are calculated for each patch. The patches are classified as follows:

  If $\lambda_{2} < c_{1}$ and $\lambda_{3} < c_{2}$ then face classified as Plane;
  If $\lambda_{2} < c_{3}$ and $\lambda_{3} < c_{4}$ then face classified as Edge;
  If $\lambda_{3} > c_{5}$ then face classified as Corner;
  In all other cases face is classified as Transitional;

  - string modelPath -- absolute input .obj path. Must be specified;
  - string meshNoisedDir -- absolute output .obj path. Must be specified;
  - double c1 -- Patch type parameter. Default value is 0.005;
  - double c2 -- Patch type parameter. Default value is 0.01;
  - double c3 -- Patch type parameter. Default value is 0.01;
  - double c4 -- Patch type parameter. Default value is 0.01;
  - double c5 -- Patch type parameter. Default value is 0.05;
  - int type -- identifies types of faces to be removed. Possible values are: 1 -- Plane, 2 -- Edge, 3 -- Corner, 4 -- Transition, 5 -- Plane and Edge, 6 -- Edge and Corner, 7 -- Plane and Corner, 8 -- Plane, Edge and Corner. Default value is 3;
  - double percentage -- defines percentage of detected faces of specified type(s) to be deleted. The interval of definition is $\[0.0, 100.0\]$. Default value is 30.0;
  - int ringsNumber -- the number of face rings in the vicinity of each face when Patch type is defined. Default value is 2;
  - int ringsNumberToDelete -- the number of face rings that will be deleted with each face that will be selected for removing. Default value is 4;
  - int regionRadius -- the number of faces theat will be considered for Patch type calculation is the vicinity of face. Default value is 16;
  - int seed -- random component determination. Default value is 0;


## Examples

```bash
import CPPNoiseGenerator
modelPath = 'absolute_path/model.obj'
```

```bash
CPPNoiseGenerator.GenerateGaussianNoisePy(modelPath, 0.2, 0)
CPPNoiseGenerator.GenerateImpulsiveNoisePy(modelPath, 0.2, 0.5, 0)
CPPNoiseGenerator.GenerateExponentialNoisePy(modelPath, 3, 0)
CPPNoiseGenerator.GenerateExtremeValueNoisePy(modelPath, 0.0, 0.3, 0)
CPPNoiseGenerator.GenerateGammaNoisePy(modelPath, 0.1, 0.2, 0)
CPPNoiseGenerator.GenerateLaplaceNoisePy(modelPath, 0.0, 0.3, 0)
CPPNoiseGenerator.GenerateLogNormalNoisePy(modelPath, -5, 40, 0)
CPPNoiseGenerator.GenerateUniformNoisePy(modelPath, 0.2, 0.4, 0)
CPPNoiseGenerator.GenerateWeibullNoisePy(modelPath, 1, 0.2, 0)
```

```bash
import CPPNoiseGenerator
modelPath = 'absolute_path/model.obj'
modelNoisedPath = 'absolute_path/model_noised.obj'
```

```bash
CPPNoiseGenerator.GenerateRandomVerticesTopologyNoisePy(modelPath, modelNoisedPath, 5)
CPPNoiseGenerator.GenerateOneRandomClusterTopologyNoisePy(modelPath, modelNoisedPath, 0.5, 50, 40.0)
CPPNoiseGenerator.GenerateSetOfRandomClustersTopologyNoisePy(modelPath, modelNoisedPath, 0.5, 50, 40.0, 4, 40, 0)
CPPNoiseGenerator.GenerateSetOfRandomClustersTopologyNoisePy(modelPath, modelNoisedPath, 0.5, 30, 20.0, 5, 40, 1)
CPPNoiseGenerator.GenerateSetPatchTopologyNoisePy(modelPath, modelNoisedPath, 0.005, 0.01, 0.01, 0.01, 0.05, 3, 30, 2, 3, 16)
```



