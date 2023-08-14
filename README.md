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
  
* `GenerateImpulsiveNoisePy(string modelPath, double sigma, double verticesPortion, int noiseDirection, int seed)`: sigma -- standard deviation of the Gaussian noise model, verticesPortion -- a portion of vertices to be noised,

- string modelPath -- absolute .obj path. Must be specified;
- double sigma -- standard deviation of the Gaussian noise model. Defaultvalue is 0.2;
- verticesPortion -- a portion of vertices to be noised. Default value is 0.2;
- int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
- int seed -- random component determination. Default value is 0;
  
* `GenerateExponentialNoisePy(string meshDir, double lambda, int noiseDirection, int seed)` -- generates node noise distributed by Exponential PDF: $P(x | \lambda) = \lambda e^{-\lambda x}$.

- string modelPath -- absolute .obj path. Must be specified;
- double lambda -- Exponential distribution $\lambda$ rate parameter. Default value is 7.0;
- int noiseDirection -- the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Default value is 0;
- int seed -- random component determination. Default value is 0;
  
* GenerateExtremeValueNoisePy(string meshDir, double a, double b, int noiseDirection, int seed)

- string modelPath -- absolute .obj path. Must be specified;
- double sigma -- standard deviation of the Gaussian noise model. Defaultvalue is 0.2;
- double lambda -- Exponential distribution parameter. Default value is 7.0;
- int seed -- random component determination. Default value is 0;
  
* ikpoikp
* dfd
* dfdf
* dfdfd
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



