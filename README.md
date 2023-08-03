# Project description

The project was created in order to provide an opportunity to generate noise for 3D meshes. This is necessary to create a dataset for learning-based denoising algorithms. This project provides a comprehensive library with 3D mesh noise generation methods.

The project is divided into two subprojects:

* `CPPNoiseGenerator_PyLib`: the noise generation algorithms are implemented in this part. All algorithms are packed into dll library.
* `NoiseGenLib`: the Python interface for the dll library is created in this part. All methods are packed into Python module.


# Algorithms provided

The algorithms implemented are divided into two node noise algorithms and topology noise algorithms. Each function reads an obj file from specified directory and writes another obj file to the directory where original obj file is placed.

Each node noise generation function has a modelPath argument of string type which specifies an absolute obj path and noiseDirection argument of int type which specifies the direction of nodes shift (0 -- node normal direction, 1 -- random direction). Node normal is calculated as a mean of adjacent faces normals.

* `GenerateGaussianNoisePy(string modelPath, double sigma, int noiseDirection)`: sigma -- dtandard deviation of the Gaussian noise model;
* `GenerateImpulsiveNoisePy(stringmodelPath, double sigma, double verticesPortion, int noiseDirection)`: 


# Examples

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



