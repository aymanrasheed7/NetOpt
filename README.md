# NetOpt
## General command
```
outFile ARG1 ARG2 ARG3 OPTIONALARG4 OPTIONALARG5
```
## Argument examples
```
ARG1: "E12" (using {1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2}), "ONE", "INT", "ODD"
ARG2: "E", "PI", "PHI", "SQRT" (sqrt(nR)), "SQRT2.5", "2.5", "-2.5" (2.5 with optional exclution)
ARG3: "1", "2", ..., "12" (nR = number of resistors)
OPTIONALARG4: "1", "2", ..., "12" (table size controller)
OPTIONALARG5: "1", "2", ..., "256" (number of threads to create)
```
## Example usage
### Windows
```bash
g++ .\NetOptThread.cpp; .\a.exe E12 SQRT 10
```
Or
```bash
g++ -lopencl .\NetOptGPU.cpp; .\a.exe E12 SQRT 10
```
### Linux
```bash
g++ -pthread NetOptThread.cpp; ./a.out E12 SQRT 10
```
Or
```bash
g++ -pthread -lopencl NetOptGPU.cpp; ./a.out E12 SQRT 10
```
## Example output
### Windows
```
Solution: (1.2+(3.9|(2.2+(((1.0|(1.5+4.7))|(2.7|5.6))+(1.8|3.3)))))
Result: 3.1622776652262714 (845658744/267420775)
Target: 3.1622776601683795
Cost: 5.057892e-09
CPU time (s): 12.819000000
Execution time (s): 12.815042700
```
### Linux
```
Solution: (1.2+(3.9|(2.2+(((1.0|(1.5+4.7))|(2.7|5.6))+(1.8|3.3)))))
Result: 3.1622776652262714 (845658744/267420775)
Target: 3.1622776601683795
Cost: 5.057892e-09
CPU time (s): 135.711018000
Execution time (s): 8.269306404
```
## About GPU
Performance on GPU was worse than that on CPU. To run on GPU proper drivers are needed.

Tested using OpenCL 3.0 CUDA driver on NVIDIA GeForce RTX 3050 Ti Laptop GPU in Windows.

Did not work when using OpenCL 3.0 NEO driver on Intel(R) Iris(R) Xe Graphics in Windows.
